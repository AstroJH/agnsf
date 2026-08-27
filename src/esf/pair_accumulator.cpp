#include <cstddef>
#include <stdexcept>
#include <vector>

#include <esf/pair_accumulator.hpp>

namespace agnsf {
namespace esf {

namespace {

void validate_redshift(double redshift)
{
    if (redshift <= -1.0) {
        throw std::invalid_argument(
            "redshift must be > -1"
        );
    }
}


std::vector<BinAccumulator> accumulate_light_curve_impl(
    const double* time,
    const double* value,
    const double* error,
    std::size_t n,
    const LagBins& bins,
    double redshift
)
{
    // TODO: Consider validating redshift at the user-facing layer.
    validate_redshift(redshift);
    const double lag_scale = 1.0 / (1.0 + redshift);

    std::vector<BinAccumulator> accumulators(
        bins.size()
    );

    const double min_lag = bins.min();
    const double max_lag = bins.max();

    // The input time array is required to be sorted.
    // Therefore, for fixed i, lag increases monotonically
    // as j increases.
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i + 1; j < n; ++j) {

            const double lag =
                (time[j] - time[i]) * lag_scale;

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

            const double delta = value[j] - value[i];

            accumulators[bin].add(
                delta,
                error[i],
                error[j]
            );
        }
    }

    return accumulators;
}

} // namespace


std::vector<BinAccumulator> accumulate_light_curve(
    const agnsf::LightCurve& light_curve,
    const LagBins& bins,
    double redshift
)
{
    return accumulate_light_curve_impl(
        light_curve.time_data(),
        light_curve.value_data(),
        light_curve.error_data(),
        light_curve.size(),
        bins,
        redshift
    );
}


std::vector<BinAccumulator> accumulate_light_curve(
    const agnsf::LightCurveView& light_curve,
    const LagBins& bins,
    double redshift
)
{
    return accumulate_light_curve_impl(
        light_curve.time_data(),
        light_curve.value_data(),
        light_curve.error_data(),
        light_curve.size(),
        bins,
        redshift
    );
}

} // namespace esf
} // namespace agnsf