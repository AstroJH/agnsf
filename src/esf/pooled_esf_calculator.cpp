#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <memory>
#include <random>
#include <stdexcept>
#include <thread>
#include <utility>
#include <vector>

#include <esf/bin_accumulator.hpp>
#include "esf_uncertainty.hpp"
#include <esf/pair_accumulator.hpp>
#include <esf/pooled_esf_calculator.hpp>
#include <esf/sf_uncertainty_estimator.hpp>

namespace agnsf {
namespace esf {
namespace
{

void validate_config(
    const UncertaintyConfig& config
)
{
    if (config.measurement != UncertaintyMethod::Off &&
        config.measurement != UncertaintyMethod::Analytic) {

        throw std::invalid_argument(
            "PooledESFCalculator: measurement uncertainty supports "
            "only Off or Analytic"
        );
    }

    if (config.sampling == UncertaintyMethod::Analytic) {

        throw std::invalid_argument(
            "PooledESFCalculator: analytic sampling uncertainty is "
            "not defined for pooled ESF; use Jackknife or Bootstrap"
        );
    }

    if (config.sampling == UncertaintyMethod::Bootstrap &&
        config.n_bootstrap == 0) {

        throw std::invalid_argument(
            "PooledESFCalculator: n_bootstrap must be >= 1"
        );
    }
}


/**
 * Fill the measurement uncertainty of every pooled bin using the
 * analytic within-bin standard-error estimator.
 */
void fill_measurement(
    std::vector<SFBinResult>& results,
    const std::vector<BinAccumulator>& accumulators,
    SFMethod method
)
{
    std::unique_ptr<SFMeasurementUncertaintyEstimator> estimator =
        make_sf_measurement_uncertainty_estimator(method);

    for (std::size_t j = 0; j < results.size(); ++j) {
        results[j].measurement =
            estimator->estimate(accumulators[j]);
    }
}


/**
 * Fill the sampling uncertainty of every pooled bin.
 *
 * The pooled ESF has no per-curve statistics, so the resampling is
 * done on the retained per-curve BinAccumulators:
 *
 *   Jackknife: leave-one-curve-out pooled SF per bin;
 *   Bootstrap: per-bin resample of curves with replacement.
 *
 * Each bin uses an independent bootstrap resample (the same seed
 * reproduces the result).
 */
void fill_sampling(
    std::vector<SFBinResult>& results,
    const std::vector<BinAccumulator>& accumulators,
    const std::vector<std::vector<BinAccumulator>>& curve_stats,
    SFMethod method,
    const UncertaintyConfig& config
)
{
    const std::size_t n_curves = curve_stats.size();

    for (std::size_t j = 0; j < results.size(); ++j) {

        if (!std::isfinite(results[j].sf)) {
            continue;
        }

        std::vector<double> values;
        values.reserve(n_curves);

        switch (config.sampling) {

            case UncertaintyMethod::Jackknife: {

                // Leave-one-curve-out pooled SF.
                for (std::size_t k = 0; k < n_curves; ++k) {

                    BinAccumulator leave_one_out =
                        accumulators[j];

                    leave_one_out.subtract(
                        curve_stats[k][j]
                    );

                    values.push_back(
                        leave_one_out.sf(method)
                    );
                }

                results[j].sampling =
                    detail::jackknife_interval(
                        values,
                        false
                    );
                break;
            }

            case UncertaintyMethod::Bootstrap: {

                std::mt19937 rng(config.bootstrap_seed);
                std::uniform_int_distribution<std::size_t>
                    index(0, n_curves - 1);

                for (std::size_t r = 0;
                     r < config.n_bootstrap;
                     ++r) {

                    BinAccumulator resampled;

                    for (std::size_t k = 0; k < n_curves; ++k) {
                        resampled.merge(
                            curve_stats[index(rng)][j]
                        );
                    }

                    values.push_back(
                        resampled.sf(method)
                    );
                }

                results[j].sampling =
                    detail::bootstrap_interval(
                        values,
                        config.n_bootstrap,
                        config.bootstrap_seed,
                        false
                    );
                break;
            }

            case UncertaintyMethod::Off:
            case UncertaintyMethod::Analytic:
                break;
        }
    }
}


/**
 * Calculate pooled ESF from light-curve views.
 *
 * All pair contributions from all light curves are accumulated
 * into the same lag bins.
 *
 * The light-curve data themselves are never copied.
 */
SFResult calculate_pooled(
    const std::vector<agnsf::LightCurveView>& data,
    const LagBins& bins,
    SFMethod method,
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

    const bool want_sampling =
        config.sampling == UncertaintyMethod::Jackknife ||
        config.sampling == UncertaintyMethod::Bootstrap;

    /*
     * Each worker processes complete light curves and owns its
     * accumulator array. Therefore, workers never write to the
     * same BinAccumulator and no synchronization is required
     * during pair accumulation.
     */

    if (data.empty()) {
        return SFResult(
            std::vector<SFBinResult>(bins.size())
        );
    }

    const std::size_t hardware_threads =
        std::thread::hardware_concurrency();

    const std::size_t num_threads =
        std::max<std::size_t>(
            1,
            std::min(
                hardware_threads == 0
                    ? std::size_t{1}
                    : hardware_threads,
                data.size()
            )
        );

    /*
     * One accumulator array per worker.
     *
     * local_accumulators[t][b]
     *     = statistics collected by worker t
     *       for lag bin b.
     */
    std::vector<std::vector<BinAccumulator>>
        local_accumulators(
            num_threads,
            std::vector<BinAccumulator>(bins.size())
        );

    /*
     * Per-curve accumulated statistics retained for the jackknife /
     * bootstrap sampling. Each curve is written by exactly one
     * worker, so no synchronization is needed.
     */
    std::vector<std::vector<BinAccumulator>> curve_stats;

    if (want_sampling) {
        curve_stats.assign(
            data.size(),
            std::vector<BinAccumulator>(bins.size())
        );
    }

    std::atomic<std::size_t> next_index{0};

    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    for (std::size_t worker_id = 0; worker_id < num_threads; ++worker_id) {

        workers.emplace_back(
            [&data,
             &bins,
             &next_index,
             &local_accumulators,
             &curve_stats,
             want_sampling,
             redshift,
             redshifts,
             worker_id]()
            {
                /*
                 * Dynamically acquire the next light curve.
                 *
                 * This avoids assigning a fixed number of light
                 * curves to each worker, which can lead to poor
                 * load balancing when light curves have very
                 * different numbers of observations.
                 */
                while (true) {

                    const std::size_t index =
                        next_index.fetch_add(
                            1,
                            std::memory_order_relaxed
                        );

                    if (index >= data.size()) {
                        break;
                    }

                    // Per-curve redshift when a vector was supplied.
                    const double curve_redshift =
                        redshifts != nullptr
                            ? (*redshifts)[index]
                            : redshift;

                    const auto result =
                        accumulate_light_curve(
                            data[index],
                            bins,
                            curve_redshift
                        );

                    if (want_sampling) {
                        curve_stats[index] = result;
                    }

                    /*
                     * Each worker writes only to its own entry,
                     * so no lock is needed here.
                     */
                    for (std::size_t bin = 0; bin < bins.size(); ++bin) {

                        local_accumulators[worker_id][bin]
                            .merge(result[bin]);
                    }
                }
            }
        );
    }

    for (auto& worker : workers) {
        worker.join();
    }

    /*
     * Merge all worker-local accumulators into one global
     * accumulator array.
     *
     * This is the pair-level pooling step:
     *
     *   global_bin = sum(worker_bin)
     *
     * over all workers.
     */
    std::vector<BinAccumulator> accumulators(
        bins.size()
    );

    for (const auto& worker_accumulators :
         local_accumulators) {

        for (std::size_t bin = 0;
             bin < bins.size();
             ++bin) {

            accumulators[bin].merge(
                worker_accumulators[bin]
            );
        }
    }

    /*
     * Convert the final accumulated pair statistics into
     * SFBinResult objects.
     */
    std::vector<SFBinResult> results;
    results.reserve(bins.size());

    for (const auto& accumulator : accumulators) {

        SFBinResult result;

        result.count = accumulator.count();
        result.sf_squared = accumulator.sf_squared(method);
        result.sf = accumulator.sf(method);

        results.push_back(result);
    }

    if (config.measurement != UncertaintyMethod::Off) {
        fill_measurement(results, accumulators, method);
    }

    if (want_sampling) {
        fill_sampling(
            results,
            accumulators,
            curve_stats,
            method,
            config
        );
    }

    return SFResult(
        std::move(results)
    );
}
} // namespace


SFResult PooledESFCalculator::calculate(
    const std::vector<agnsf::LightCurve>& data,
    const LagBins& bins,
    SFMethod method,
    const UncertaintyConfig& config,
    double redshift
) const
{
    /*
     * Convert owning light curves into non-owning views.
     */
    std::vector<agnsf::LightCurveView> views;
    views.reserve(data.size());

    for (const auto& light_curve : data) {
        views.push_back(light_curve.view());
    }

    return calculate_pooled(
        views,
        bins,
        method,
        config,
        redshift,
        nullptr
    );
}


SFResult PooledESFCalculator::calculate(
    const std::vector<agnsf::LightCurveView>& data,
    const LagBins& bins,
    SFMethod method,
    const UncertaintyConfig& config,
    double redshift
) const
{
    return calculate_pooled(
        data,
        bins,
        method,
        config,
        redshift,
        nullptr
    );
}


SFResult PooledESFCalculator::calculate(
    const std::vector<agnsf::LightCurve>& data,
    const LagBins& bins,
    SFMethod method,
    const UncertaintyConfig& config,
    const std::vector<double>& redshifts
) const
{
    std::vector<agnsf::LightCurveView> views;
    views.reserve(data.size());

    for (const auto& light_curve : data) {
        views.push_back(light_curve.view());
    }

    return calculate_pooled(
        views,
        bins,
        method,
        config,
        0.0,
        &redshifts
    );
}


SFResult PooledESFCalculator::calculate(
    const std::vector<agnsf::LightCurveView>& data,
    const LagBins& bins,
    SFMethod method,
    const UncertaintyConfig& config,
    const std::vector<double>& redshifts
) const
{
    return calculate_pooled(
        data,
        bins,
        method,
        config,
        0.0,
        &redshifts
    );
}

} // namespace esf
} // namespace agnsf
