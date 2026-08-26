#include "esf_uncertainty.hpp"  // IWYU pragma: private

#include <algorithm>
#include <cmath>
#include <random>

namespace agnsf {
namespace esf {
namespace detail {

namespace {

double sum_of(const std::vector<double>& values)
{
    double sum = 0.0;

    for (const double value : values) {
        sum += value;
    }

    return sum;
}


double mean_of(const std::vector<double>& values)
{
    return sum_of(values) / static_cast<double>(values.size());
}


/**
 * Apply the per-replicate statistic: identity (MeanSf) or
 * sqrt(mean) with the argument clamped at 0 (SqrtMeanSquared, so
 * that a noise-dominated bin still yields a real sf value).
 */
double replicate_statistic(
    double mean,
    bool sqrt_mean
)
{
    if (!sqrt_mean) {
        return mean;
    }

    return std::sqrt(std::max(mean, 0.0));
}

} // namespace


SFUncertainty analytic_interval(
    const std::vector<double>& values
)
{
    SFUncertainty result;

    if (values.size() < 2) {
        return result;
    }

    const double n =
        static_cast<double>(values.size());

    double sum = 0.0;
    double sum_squared = 0.0;

    for (const double value : values) {
        sum += value;
        sum_squared += value * value;
    }

    const double mean = sum / n;

    // Sample variance of the per-curve values; clamp rounding noise.
    double sample_var =
        (sum_squared - n * mean * mean) /
        (n - 1.0);

    sample_var = std::max(sample_var, 0.0);

    const double se =
        std::sqrt(sample_var / n);

    result.lower = mean - se;
    result.upper = mean + se;

    return result;
}


SFUncertainty jackknife_interval(
    const std::vector<double>& values,
    bool sqrt_mean
)
{
    SFUncertainty result;

    const std::size_t n_size = values.size();

    if (n_size < 2) {
        return result;
    }

    const double n =
        static_cast<double>(n_size);

    const double total = sum_of(values);

    // Leave-one-out statistics theta_{(-k)}.
    std::vector<double> leave_one_out;
    leave_one_out.reserve(n_size);

    for (const double value : values) {
        const double mean_without_k =
            (total - value) / (n - 1.0);

        leave_one_out.push_back(
            replicate_statistic(mean_without_k, sqrt_mean)
        );
    }

    // Jackknife variance:
    //   var = (n - 1) / n * sum_k (theta_{(-k)} - theta_{(.)})^2
    const double leave_one_out_mean =
        mean_of(leave_one_out);

    double sum_squared_deviation = 0.0;

    for (const double value : leave_one_out) {
        const double deviation =
            value - leave_one_out_mean;

        sum_squared_deviation +=
            deviation * deviation;
    }

    const double variance =
        ((n - 1.0) / n) * sum_squared_deviation;

    const double sigma =
        std::sqrt(variance);

    // Center on the full-sample statistic.
    const double full_mean = total / n;

    const double point =
        replicate_statistic(full_mean, sqrt_mean);

    result.lower = point - sigma;
    result.upper = point + sigma;

    return result;
}


SFUncertainty bootstrap_interval(
    const std::vector<double>& values,
    std::size_t n_bootstrap,
    std::uint32_t seed,
    bool sqrt_mean
)
{
    SFUncertainty result;

    const std::size_t n = values.size();

    if (n < 2 || n_bootstrap == 0) {
        return result;
    }

    std::mt19937 rng(seed);
    std::uniform_int_distribution<std::size_t>
        index(0, n - 1);

    std::vector<double> replicates;
    replicates.reserve(n_bootstrap);

    for (std::size_t r = 0; r < n_bootstrap; ++r) {

        double sum = 0.0;

        // Resample n curves with replacement.
        for (std::size_t k = 0; k < n; ++k) {
            sum += values[index(rng)];
        }

        const double statistic =
            replicate_statistic(
                sum / static_cast<double>(n),
                sqrt_mean
            );

        if (std::isfinite(statistic)) {
            replicates.push_back(statistic);
        }
    }

    if (replicates.size() < 2) {
        return result;
    }

    std::sort(replicates.begin(), replicates.end());

    // 16th / 84th percentiles (nearest rank).
    const std::size_t lower_index =
        static_cast<std::size_t>(
            0.16 * static_cast<double>(replicates.size() - 1)
        );

    const std::size_t upper_index =
        static_cast<std::size_t>(
            0.84 * static_cast<double>(replicates.size() - 1)
        );

    result.lower = replicates[lower_index];
    result.upper = replicates[upper_index];

    return result;
}


SFUncertainty map_interval_to_sf(
    const SFUncertainty& interval,
    bool sqrt_transform
)
{
    if (!sqrt_transform) {
        return interval;
    }

    if (!interval.estimated()) {
        return SFUncertainty{};
    }

    SFUncertainty result;

    result.lower =
        std::sqrt(std::max(interval.lower, 0.0));

    result.upper =
        std::sqrt(std::max(interval.upper, 0.0));

    return result;
}

} // namespace detail
} // namespace esf
} // namespace agnsf
