#include <cstddef>
#include <vector>

#include <esf/pair_accumulator.hpp>

namespace esf {

std::vector<BinAccumulator> accumulate_light_curve(
    const LightCurve& light_curve,
    const LagBins& bins
)
{
    std::vector<BinAccumulator> accumulators(
        bins.size()
    );

    const std::size_t n =
        light_curve.size();

    const double* time =
        light_curve.time_data();

    const double* value =
        light_curve.value_data();

    const double* error =
        light_curve.error_data();

    const double min_lag =
        bins.min();

    const double max_lag =
        bins.max();

    // The input time array is required to be sorted.
    // Therefore, for fixed i, lag increases monotonically
    // as j increases.
    for (std::size_t i = 0; i < n; ++i) {

        for (std::size_t j = i + 1; j < n; ++j) {

            const double lag =
                time[j] - time[i];

            // Since time is sorted, all subsequent j values
            // have equal or larger lag.
            if (lag >= max_lag) {
                break;
            }

            if (lag < min_lag) {
                continue;
            }

            std::size_t bin;

            if (!bins.try_index(lag, bin)) {
                continue;
            }

            const double delta =
                value[j] - value[i];

            accumulators[bin].add(
                delta,
                error[i],
                error[j]
            );
        }
    }

    return accumulators;
}

}