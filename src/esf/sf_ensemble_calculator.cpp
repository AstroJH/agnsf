#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

#include "esf_uncertainty.hpp"
#include <esf/sf_calculator.hpp>
#include <esf/sf_ensemble_calculator.hpp>

namespace agnsf {
namespace esf {

namespace {

/**
 * Turn accumulated per-bin statistics into an SFResult (point
 * estimate only; uncertainty fields are filled separately).
 *
 * For SqrtMeanSquared, sum holds the sum of per-curve SF^2 values
 * and the result is:
 *
 *   SF^2 = sum / count
 *   SF   = sqrt(SF^2)
 *
 * For MeanSf, sum holds the sum of per-curve SF values and the
 * result is:
 *
 *   SF   = sum / count
 *   SF^2 = SF^2
 */
std::vector<SFBinResult> make_result(
    const std::vector<double>& sum,
    const std::vector<std::size_t>& count,
    SFEnsembleCalculator::Method method
)
{
    std::vector<SFBinResult> results;
    results.reserve(sum.size());

    for (std::size_t i = 0;
         i < sum.size();
         ++i) {

        SFBinResult result;
        result.count = count[i];

        if (count[i] == 0) {

            result.sf_squared =
                std::numeric_limits<double>::quiet_NaN();

            result.sf =
                std::numeric_limits<double>::quiet_NaN();

        } else if (method == SFEnsembleCalculator::Method::MeanSf) {

            result.sf =
                sum[i] /
                static_cast<double>(count[i]);

            result.sf_squared =
                result.sf * result.sf;

        } else {

            result.sf_squared =
                sum[i] /
                static_cast<double>(count[i]);

            result.sf =
                std::sqrt(result.sf_squared);
        }

        results.push_back(result);
    }

    return results;
}


/**
 * True when the aggregated working scale is SF^2 (SqrtMeanSquared);
 * for MeanSf the working scale is SF itself.
 */
bool uses_squared_scale(
    SFEnsembleCalculator::Method method
)
{
    return method ==
        SFEnsembleCalculator::Method::SqrtMeanSquared;
}


/**
 * Working-scale sigma of one curve's interval half-width:
 *
 *   MeanSf:          sigma = half_width
 *   SqrtMeanSquared: sigma(sf^2) ~= 2 * sf * half_width  (delta method)
 */
double working_sigma_of(
    double half_width,
    double sf,
    SFEnsembleCalculator::Method method
)
{
    return uses_squared_scale(method)
        ? 2.0 * sf * half_width
        : half_width;
}


void validate_config(
    const UncertaintyConfig& config
)
{
    if (config.measurement != UncertaintyMethod::Off &&
        config.measurement != UncertaintyMethod::Analytic &&
        config.measurement != UncertaintyMethod::MonteCarlo) {

        throw std::invalid_argument(
            "SFEnsembleCalculator: measurement uncertainty supports "
            "Off, Analytic, or MonteCarlo"
        );
    }

    if (config.within != UncertaintyMethod::Off &&
        config.within != UncertaintyMethod::Analytic) {

        throw std::invalid_argument(
            "SFEnsembleCalculator: within uncertainty supports "
            "Off or Analytic"
        );
    }

    if (config.sampling == UncertaintyMethod::Bootstrap &&
        config.n_bootstrap == 0) {

        throw std::invalid_argument(
            "SFEnsembleCalculator: n_bootstrap must be >= 1"
        );
    }
}


/**
 * Combine per-curve half-widths in quadrature (independent curves)
 * and map the resulting interval on the working scale (SF or SF^2)
 * onto sf.
 */
std::vector<Uncertainty> combine_working_intervals(
    const std::vector<double>& sum2,
    const std::vector<std::size_t>& count,
    const std::vector<SFBinResult>& results,
    SFEnsembleCalculator::Method method
)
{
    const bool sqrt_transform =
        uses_squared_scale(method);

    std::vector<Uncertainty> intervals(
        results.size()
    );

    for (std::size_t j = 0; j < results.size(); ++j) {

        const double n =
            static_cast<double>(count[j]);

        if (n < 1.0 ||
            !std::isfinite(results[j].sf)) {
            continue;
        }

        const double sigma =
            std::sqrt(sum2[j]) / n;

        // Working-scale mean (SF for MeanSf, SF^2 for SqrtMeanSquared).
        const double working_mean =
            sqrt_transform
                ? results[j].sf_squared
                : results[j].sf;

        Uncertainty working;
        working.lower = working_mean - sigma;
        working.upper = working_mean + sigma;

        intervals[j] =
            detail::map_interval_to_sf(
                working,
                sqrt_transform
            );
    }

    return intervals;
}

/**
 * Fill the measurement uncertainty of every bin.
 *
 * Per-curve measurement uncertainties (half-widths on sf) are
 * combined in quadrature under the independence assumption:
 *
 *   MeanSf:            sigma^2 = sum(hw_k^2) / n_meas^2
 *   SqrtMeanSquared:   sigma(sf^2)_k ~= 2 * sf_k * hw_k   (delta method)
 *
 * The combined sigma is applied on the working scale and mapped to sf
 * (sqrt mapping makes the result asymmetric for SqrtMeanSquared).
 */
void fill_measurement(
    std::vector<SFBinResult>& results,
    const std::vector<double>& meas2,
    const std::vector<std::size_t>& meas_count,
    SFEnsembleCalculator::Method method
)
{
    const std::vector<Uncertainty> intervals =
        combine_working_intervals(
            meas2,
            meas_count,
            results,
            method
        );

    for (std::size_t j = 0; j < results.size(); ++j) {
        results[j].measurement = intervals[j];
    }
}


void fill_within(
    std::vector<SFBinResult>& results,
    const std::vector<double>& within2,
    const std::vector<std::size_t>& within_count,
    SFEnsembleCalculator::Method method
)
{
    const std::vector<Uncertainty> intervals =
        combine_working_intervals(
            within2,
            within_count,
            results,
            method
        );

    for (std::size_t j = 0; j < results.size(); ++j) {
        results[j].within = intervals[j];
    }
}


/**
 * Fill the sampling uncertainty of every bin.
 *
 * All methods operate on the per-curve values (SF for MeanSf, SF^2
 * for SqrtMeanSquared) and produce an interval on sf:
 *
 *   Analytic:  mean +/- std/sqrt(n) on the working scale, mapped to sf
 *   Jackknife: leave-one-curve-out, sigma centered on the point value
 *   Bootstrap: 16/84 percentile interval of resampled statistics
 */
void fill_sampling(
    std::vector<SFBinResult>& results,
    const std::vector<double>& sum,
    const std::vector<double>& sum2,
    const std::vector<std::size_t>& count,
    const std::vector<std::vector<double>>& values,
    SFEnsembleCalculator::Method method,
    const UncertaintyConfig& config
)
{
    const bool sqrt_transform =
        uses_squared_scale(method);

    for (std::size_t j = 0; j < results.size(); ++j) {

        if (!std::isfinite(results[j].sf)) {
            continue;
        }

        Uncertainty interval;

        switch (config.sampling) {

            case UncertaintyMethod::Analytic: {

                const double n =
                    static_cast<double>(count[j]);

                if (count[j] < 2) {
                    break;
                }

                const double mean = sum[j] / n;

                double sample_var =
                    (sum2[j] - n * mean * mean) /
                    (n - 1.0);

                sample_var =
                    std::max(sample_var, 0.0);

                const double se =
                    std::sqrt(sample_var / n);

                interval.lower = mean - se;
                interval.upper = mean + se;

                interval =
                    detail::map_interval_to_sf(
                        interval,
                        sqrt_transform
                    );
                break;
            }

            case UncertaintyMethod::Jackknife: {

                // Collect the finite per-curve values of this bin.
                std::vector<double> bin_values;

                for (const auto& row : values) {
                    const double value = row[j];

                    if (std::isfinite(value)) {
                        bin_values.push_back(value);
                    }
                }

                interval =
                    detail::jackknife_interval(
                        bin_values,
                        sqrt_transform
                    );
                break;
            }

            case UncertaintyMethod::Bootstrap: {

                std::vector<double> bin_values;

                for (const auto& row : values) {
                    const double value = row[j];

                    if (std::isfinite(value)) {
                        bin_values.push_back(value);
                    }
                }

                interval =
                    detail::bootstrap_interval(
                        bin_values,
                        config.n_bootstrap,
                        config.bootstrap_seed,
                        sqrt_transform
                    );
                break;
            }

            case UncertaintyMethod::Off:
            case UncertaintyMethod::MonteCarlo:
                break;
        }

        results[j].sampling = interval;
    }
}


SFResult calculate_ensemble(
    const std::vector<agnsf::LightCurveView>& data,
    const LagBins& bins,
    SFMethod sf_method,
    SFEnsembleCalculator::Method method,
    const UncertaintyConfig& config,
    double redshift,
    const std::vector<double>* redshifts  // nullptr => scalar redshift
)
{
    validate_config(config);

    if (redshifts != nullptr &&
        redshifts->size() != data.size()) {

        throw std::invalid_argument(
            "redshifts size must match the number of light curves"
        );
    }

    if (redshift <= -1.0) {
        throw std::invalid_argument(
            "redshift must be > -1"
        );
    }

    const bool want_measurement =
        config.measurement != UncertaintyMethod::Off;

    const bool want_within =
        config.within == UncertaintyMethod::Analytic;

    const bool want_analytic_sampling =
        config.sampling == UncertaintyMethod::Analytic;

    const bool want_resampling =
        config.sampling == UncertaintyMethod::Jackknife ||
        config.sampling == UncertaintyMethod::Bootstrap;

    if (data.empty()) {

        const std::vector<double> sum(
            bins.size(),
            0.0
        );

        const std::vector<std::size_t> count(
            bins.size(),
            0
        );

        return SFResult(
            make_result(sum, count, method)
        );
    }


    std::size_t n_threads =
        std::thread::hardware_concurrency();

    if (n_threads == 0) {
        n_threads = 1;
    }

    n_threads =
        std::min(n_threads, data.size());


    /*
     * Each thread has its own partial results.
     *
     * thread_sum / thread_count      : point estimate (existing)
     * thread_sum2                    : sum of squared per-curve values
     *                                  (analytic sampling)
     * thread_meas2 / thread_meas_count : per-curve measurement
     *                                  variance sum (measurement)
     */
    std::vector<std::vector<double>> thread_sum(
        n_threads,
        std::vector<double>(bins.size(), 0.0)
    );

    std::vector<std::vector<std::size_t>> thread_count(
        n_threads,
        std::vector<std::size_t>(bins.size(), 0)
    );

    std::vector<std::vector<double>> thread_sum2;

    if (want_analytic_sampling) {
        thread_sum2.assign(
            n_threads,
            std::vector<double>(bins.size(), 0.0)
        );
    }

    std::vector<std::vector<double>> thread_meas2;

    std::vector<std::vector<std::size_t>> thread_meas_count;

    if (want_measurement) {
        thread_meas2.assign(
            n_threads,
            std::vector<double>(bins.size(), 0.0)
        );

        thread_meas_count.assign(
            n_threads,
            std::vector<std::size_t>(bins.size(), 0)
        );
    }

    std::vector<std::vector<double>> thread_within2;

    std::vector<std::vector<std::size_t>> thread_within_count;

    if (want_within) {
        thread_within2.assign(
            n_threads,
            std::vector<double>(bins.size(), 0.0)
        );

        thread_within_count.assign(
            n_threads,
            std::vector<std::size_t>(bins.size(), 0)
        );
    }

    /*
     * Per-curve values (SF or SF^2) retained for jackknife /
     * bootstrap sampling. Each thread writes only its own rows.
     */
    std::vector<std::vector<double>> values;

    if (want_resampling) {
        values.assign(
            data.size(),
            std::vector<double>(
                bins.size(),
                std::numeric_limits<double>::quiet_NaN()
            )
        );
    }


    const std::size_t chunk_size =
        (data.size() + n_threads - 1) /
        n_threads;


    std::vector<std::thread> threads;
    threads.reserve(n_threads);


    for (std::size_t thread_id = 0;
         thread_id < n_threads;
         ++thread_id) {

        const std::size_t begin =
            thread_id * chunk_size;

        const std::size_t end =
            std::min(
                begin + chunk_size,
                data.size()
            );

        threads.emplace_back(
            [&, thread_id, begin, end]() {

                SFCalculator sf_calculator;

                // Per-curve SF keeps the measurement component;
                // sampling is handled at the ensemble level.
                const UncertaintyConfig per_curve_config{
                    config.measurement,
                    config.within,
                    UncertaintyMethod::Off,
                    0,
                    0
                };

                for (std::size_t i = begin;
                     i < end;
                     ++i) {

                    const double curve_redshift =
                        redshifts != nullptr
                            ? (*redshifts)[i]
                            : redshift;

                    const SFResult result =
                        sf_calculator.calculate(
                            data[i],
                            bins,
                            sf_method,
                            per_curve_config,
                            curve_redshift
                        );

                    for (std::size_t j = 0;
                         j < bins.size();
                         ++j) {

                        /*
                         * SqrtMeanSquared averages per-curve SF^2.
                         * MeanSf averages per-curve SF instead.
                         */
                        const double contribution =
                            method ==
                                SFEnsembleCalculator::Method::MeanSf
                                ? result.bin(j).sf
                                : result.bin(j).sf_squared;

                        if (!std::isfinite(contribution)) {
                            continue;
                        }

                        thread_sum[thread_id][j] +=
                            contribution;

                        ++thread_count[thread_id][j];

                        if (want_analytic_sampling) {
                            thread_sum2[thread_id][j] +=
                                contribution * contribution;
                        }

                        if (want_resampling) {
                            values[i][j] = contribution;
                        }

                        if (want_measurement) {

                            const Uncertainty& interval =
                                result.bin(j).measurement;

                            if (interval.estimated()) {

                                const double half_width =
                                    (
                                        interval.upper -
                                        interval.lower
                                    ) / 2.0;

                                if (std::isfinite(half_width)) {

                                    const double working_sigma =
                                        working_sigma_of(
                                            half_width,
                                            result.bin(j).sf,
                                            method
                                        );

                                    thread_meas2[thread_id][j] +=
                                        working_sigma * working_sigma;

                                    ++thread_meas_count[thread_id][j];
                                }
                            }
                        }

                        if (want_within) {

                            const Uncertainty& interval =
                                result.bin(j).within;

                            if (interval.estimated()) {

                                const double half_width =
                                    (
                                        interval.upper -
                                        interval.lower
                                    ) / 2.0;

                                if (std::isfinite(half_width)) {

                                    const double working_sigma =
                                        working_sigma_of(
                                            half_width,
                                            result.bin(j).sf,
                                            method
                                        );

                                    thread_within2[thread_id][j] +=
                                        working_sigma * working_sigma;

                                    ++thread_within_count[thread_id][j];
                                }
                            }
                        }
                    }
                }
            }
        );
    }


    for (auto& thread : threads) {
        thread.join();
    }


    /*
     * Merge thread-local results.
     */
    std::vector<double> sum(
        bins.size(),
        0.0
    );

    std::vector<std::size_t> count(
        bins.size(),
        0
    );

    std::vector<double> sum2(
        bins.size(),
        0.0
    );

    std::vector<double> meas2(
        bins.size(),
        0.0
    );

    std::vector<std::size_t> meas_count(
        bins.size(),
        0
    );

    std::vector<double> within2(
        bins.size(),
        0.0
    );

    std::vector<std::size_t> within_count(
        bins.size(),
        0
    );

    for (std::size_t thread_id = 0;
         thread_id < n_threads;
         ++thread_id) {

        for (std::size_t i = 0;
             i < bins.size();
             ++i) {

            sum[i] +=
                thread_sum[thread_id][i];

            count[i] +=
                thread_count[thread_id][i];

            if (want_analytic_sampling) {
                sum2[i] +=
                    thread_sum2[thread_id][i];
            }

            if (want_measurement) {
                meas2[i] +=
                    thread_meas2[thread_id][i];

                meas_count[i] +=
                    thread_meas_count[thread_id][i];
            }

            if (want_within) {
                within2[i] +=
                    thread_within2[thread_id][i];

                within_count[i] +=
                    thread_within_count[thread_id][i];
            }
        }
    }


    // Fill uncertainty fields on the result bins.
    std::vector<SFBinResult> bins_out =
        make_result(sum, count, method);

    if (want_measurement) {
        fill_measurement(
            bins_out,
            meas2,
            meas_count,
            method
        );
    }

    if (want_within) {
        fill_within(
            bins_out,
            within2,
            within_count,
            method
        );
    }

    if (config.sampling != UncertaintyMethod::Off) {
        fill_sampling(
            bins_out,
            sum,
            sum2,
            count,
            values,
            method,
            config
        );
    }

    return SFResult(std::move(bins_out));
}

} // namespace


SFResult SFEnsembleCalculator::calculate(
    const std::vector<agnsf::LightCurve>& data,
    const LagBins& bins,
    SFMethod sf_method,
    Method method,
    const UncertaintyConfig& config,
    double redshift
) const
{
    std::vector<agnsf::LightCurveView> views;
    views.reserve(data.size());

    for (const auto& light_curve : data) {
        views.push_back(light_curve.view());
    }

    return calculate(
        views,
        bins,
        sf_method,
        method,
        config,
        redshift
    );
}


SFResult SFEnsembleCalculator::calculate(
    const std::vector<agnsf::LightCurveView>& data,
    const LagBins& bins,
    SFMethod sf_method,
    Method method,
    const UncertaintyConfig& config,
    double redshift
) const
{
    return calculate_ensemble(
        data,
        bins,
        sf_method,
        method,
        config,
        redshift,
        nullptr
    );
}


SFResult SFEnsembleCalculator::calculate(
    const std::vector<agnsf::LightCurve>& data,
    const LagBins& bins,
    SFMethod sf_method,
    Method method,
    const UncertaintyConfig& config,
    const std::vector<double>& redshifts
) const
{
    std::vector<agnsf::LightCurveView> views;
    views.reserve(data.size());

    for (const auto& light_curve : data) {
        views.push_back(light_curve.view());
    }

    return calculate(
        views,
        bins,
        sf_method,
        method,
        config,
        redshifts
    );
}


SFResult SFEnsembleCalculator::calculate(
    const std::vector<agnsf::LightCurveView>& data,
    const LagBins& bins,
    SFMethod sf_method,
    Method method,
    const UncertaintyConfig& config,
    const std::vector<double>& redshifts
) const
{
    return calculate_ensemble(
        data,
        bins,
        sf_method,
        method,
        config,
        0.0,
        &redshifts
    );
}

} // namespace esf
} // namespace agnsf
