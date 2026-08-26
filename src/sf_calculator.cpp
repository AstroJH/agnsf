#include <esf/sf_calculator.hpp>

#include <esf/bin_accumulator.hpp>

namespace esf {

namespace {

SFResult calculate_impl(
    const double* time,
    const double* value,
    const double* error,
    std::size_t size,
    const LagBins& bins,
    SFMethod method
)
{
    std::vector<BinAccumulator> accumulators(
        bins.size()
    );

    const double min_lag = bins.min();
    const double max_lag = bins.max();

    for (std::size_t i = 0; i < size; ++i) {
        for (std::size_t j = i + 1; j < size; ++j) {

            const double lag =
                time[j] - time[i];

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

    std::vector<SFBinResult> results;
    results.reserve(bins.size());

    for (const auto& accumulator : accumulators) {

        SFBinResult result;

        result.count = accumulator.count();
        result.sf_squared = accumulator.sf_squared(method);
        result.sf = accumulator.sf(method);

        results.push_back(result);
    }

    return SFResult(
        std::move(results)
    );
}

} // namespace


SFResult SFCalculator::calculate(
    const LightCurve& data,
    const LagBins& bins,
    SFMethod method
) const
{
    return calculate_impl(
        data.time_data(),
        data.value_data(),
        data.error_data(),
        data.size(),
        bins,
        method
    );
}


SFResult SFCalculator::calculate(
    const LightCurveView& data,
    const LagBins& bins,
    SFMethod method
) const
{
    return calculate_impl(
        data.time_data(),
        data.value_data(),
        data.error_data(),
        data.size(),
        bins,
        method
    );
}

} // namespace esf
