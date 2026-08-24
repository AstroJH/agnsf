#include <utility>
#include <vector>

#include <esf/pair_accumulator.hpp>
#include <esf/sf_calculator.hpp>

namespace esf {

SFResult SFCalculator::calculate(
    const LightCurve& data,
    const LagBins& bins
) const
{
    const auto accumulators =
        accumulate_light_curve(data, bins);

    std::vector<SFBinResult> results;
    results.reserve(bins.size());

    for (const auto& accumulator : accumulators) {

        SFBinResult result;

        result.count =
            accumulator.count();

        result.sf_squared =
            accumulator.sf_squared();

        result.sf =
            accumulator.sf();

        results.push_back(result);
    }

    return SFResult(
        std::move(results)
    );
}

} // namespace esf
