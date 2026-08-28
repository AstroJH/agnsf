#include <esf/sf_calculator.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <random>
#include <stdexcept>
#include <vector>

#include <esf/bin_accumulator.hpp>
#include <esf/pair_accumulator.hpp>
#include <esf/sf_uncertainty_estimator.hpp>
#include "esf_uncertainty.hpp"

namespace agnsf {
namespace esf {

namespace {

/**
 * Validate the uncertainty configuration for a single SF.
 *
 * measurement: Analytic (sigma_i propagation, independent-pairs
 *              approximation) or MonteCarlo (observation-level
 *              perturbation).
 * within:      Analytic (naive s_X / sqrt(N_pair)).
 * sampling:    ESF-only; must be Off here.
 */
void validate_config(
    const UncertaintyConfig& config
)
{
    const auto valid_measurement =
        config.measurement == UncertaintyMethod::Off ||
        config.measurement == UncertaintyMethod::Analytic ||
        config.measurement == UncertaintyMethod::MonteCarlo;

    const auto valid_within =
        config.within == UncertaintyMethod::Off ||
        config.within == UncertaintyMethod::Analytic;

    if (!valid_measurement) {
        throw std::invalid_argument(
            "SFCalculator: measurement uncertainty supports "
            "Off, Analytic, or MonteCarlo"
        );
    }

    if (!valid_within) {
        throw std::invalid_argument(
            "SFCalculator: within uncertainty supports Off or Analytic"
        );
    }

    if (config.sampling != UncertaintyMethod::Off) {
        throw std::invalid_argument(
            "SFCalculator: sampling uncertainty is not defined "
            "for a single light curve"
        );
    }
}


SFResult calculate_impl(
    const double* time,
    const double* value,
    const double* error,
    std::size_t size,
    const LagBins& bins,
    SFMethod method,
    const UncertaintyConfig& config,
    double redshift
)
{
    validate_config(config);

    const agnsf::LightCurveView view(
        time,
        value,
        error,
        size
    );

    // Point estimate + per-bin pair statistics (rest-frame lags are
    // handled inside accumulate_light_curve).
    const std::vector<BinAccumulator> accumulators =
        accumulate_light_curve(view, bins, redshift);

    std::vector<SFBinResult> results;
    results.reserve(bins.size());

    for (const auto& accumulator : accumulators) {

        SFBinResult result;

        result.count = accumulator.count();
        result.sf_squared = accumulator.sf_squared(method);
        result.sf = accumulator.sf(method);

        results.push_back(result);
    }

    // ---- measurement uncertainty: sigma_i propagation ----
    if (config.measurement == UncertaintyMethod::Analytic) {

        std::unique_ptr<SFMeasurementUncertaintyEstimator> estimator =
            make_sf_measurement_uncertainty_estimator(method);

        if (estimator) {
            for (std::size_t j = 0; j < bins.size(); ++j) {
                results[j].measurement =
                    estimator->estimate(accumulators[j]);
            }
        }

    } else if (config.measurement == UncertaintyMethod::MonteCarlo) {

        // Observation-level perturbation: f_i* = f_i + e_i with
        // e_i ~ N(0, sigma_i^2), then the SF is recomputed per
        // realization. This naturally preserves the covariance that
        // shared observations induce between pairs. Temporary value
        // arrays are created per realization (inherent to Monte Carlo).
        const std::size_t n_realizations =
            std::max<std::size_t>(config.n_bootstrap, 1);

        std::mt19937 rng(config.bootstrap_seed);
        std::normal_distribution<double> normal(0.0, 1.0);

        std::vector<double> perturbed(size);

        std::vector<std::vector<double>> per_bin_sf(
            bins.size()
        );

        for (std::size_t r = 0; r < n_realizations; ++r) {

            for (std::size_t i = 0; i < size; ++i) {
                perturbed[i] =
                    value[i] + error[i] * normal(rng);
            }

            const agnsf::LightCurveView perturbed_view(
                time,
                perturbed.data(),
                error,
                size
            );

            const std::vector<BinAccumulator> acc =
                accumulate_light_curve(
                    perturbed_view,
                    bins,
                    redshift
                );

            for (std::size_t j = 0; j < bins.size(); ++j) {
                per_bin_sf[j].push_back(acc[j].sf(method));
            }
        }

        for (std::size_t j = 0; j < bins.size(); ++j) {
            results[j].measurement =
                detail::percentile_interval(
                    per_bin_sf[j],
                    0.16,
                    0.84
                );
        }
    }

    // ---- naive within-bin statistical uncertainty ----
    if (config.within == UncertaintyMethod::Analytic) {

        std::unique_ptr<SFWithinUncertaintyEstimator> estimator =
            make_sf_within_uncertainty_estimator(method);

        if (estimator) {
            for (std::size_t j = 0; j < bins.size(); ++j) {
                results[j].within =
                    estimator->estimate(accumulators[j]);
            }
        }
    }

    return SFResult(std::move(results));
}

} // namespace


SFResult SFCalculator::calculate(
    const agnsf::LightCurve& data,
    const LagBins& bins,
    SFMethod method,
    const UncertaintyConfig& config,
    double redshift
) const
{
    return calculate_impl(
        data.time_data(),
        data.value_data(),
        data.error_data(),
        data.size(),
        bins,
        method,
        config,
        redshift
    );
}


SFResult SFCalculator::calculate(
    const agnsf::LightCurveView& data,
    const LagBins& bins,
    SFMethod method,
    const UncertaintyConfig& config,
    double redshift
) const
{
    return calculate_impl(
        data.time_data(),
        data.value_data(),
        data.error_data(),
        data.size(),
        bins,
        method,
        config,
        redshift
    );
}

} // namespace esf
} // namespace agnsf
