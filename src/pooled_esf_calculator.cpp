#include <algorithm>
#include <atomic>
#include <cstddef>
#include <thread>
#include <utility>
#include <vector>

#include <esf/pooled_esf_calculator.hpp>
#include <esf/bin_accumulator.hpp>
#include <esf/pair_accumulator.hpp>

namespace esf {
namespace
{

/**
 * Calculate pooled ESF from light-curve views.
 *
 * All pair contributions from all light curves are accumulated
 * into the same lag bins.
 *
 * The light-curve data themselves are never copied.
 */
SFResult calculate_pooled(
    const std::vector<LightCurveView>& data,
    const LagBins& bins
)
{
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

    std::atomic<std::size_t> next_index{0};

    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    for (std::size_t worker_id = 0; worker_id < num_threads; ++worker_id) {

        workers.emplace_back(
            [&data,
             &bins,
             &next_index,
             &local_accumulators,
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

                    const auto result =
                        accumulate_light_curve(
                            data[index],
                            bins
                        );

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
        result.sf_squared = accumulator.sf_squared();
        result.sf = accumulator.sf();

        results.push_back(result);
    }

    return SFResult(
        std::move(results)
    );
}
} // namespace


SFResult PooledESFCalculator::calculate(
    const std::vector<LightCurve>& data,
    const LagBins& bins
) const
{
    /*
     * Convert owning light curves into non-owning views.
     *
     * No time/value/error data are copied.
     */
    std::vector<LightCurveView> views;
    views.reserve(data.size());

    for (const auto& light_curve : data) {
        views.emplace_back(
            light_curve.time_data(),
            light_curve.value_data(),
            light_curve.error_data(),
            light_curve.size()
        );
    }

    return calculate_pooled(
        views,
        bins
    );
}


SFResult PooledESFCalculator::calculate(
    const std::vector<LightCurveView>& data,
    const LagBins& bins
) const
{
    return calculate_pooled(
        data,
        bins
    );
}

} // namespace esf