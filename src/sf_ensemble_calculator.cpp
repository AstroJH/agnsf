#include <esf/sf_ensemble_calculator.hpp>

#include <cmath>
#include <limits>
#include <vector>

#include <esf/sf_calculator.hpp>

namespace esf {

SFResult SFEnsembleCalculator::calculate(
    const std::vector<LightCurve>& data,
    const LagBins& bins
) const
{
    SFCalculator sf_calculator;

    std::vector<double> sum_sf(
        bins.size(),
        0.0
    );

    std::vector<std::size_t> count(
        bins.size(),
        0
    );

    for (const auto& light_curve : data) {

        const SFResult result =
            sf_calculator.calculate(
                light_curve,
                bins
            );

        for (std::size_t i = 0;
             i < bins.size();
             ++i) {

            const double sf =
                result.bin(i).sf;

            if (!std::isfinite(sf)) {
                continue;
            }

            sum_sf[i] += sf;
            ++count[i];
        }
    }

    std::vector<SFBinResult> results;
    results.reserve(bins.size());

    for (std::size_t i = 0;
         i < bins.size();
         ++i) {

        SFBinResult result;

        result.count = count[i];

        if (count[i] == 0) {

            result.sf =
                std::numeric_limits<double>::quiet_NaN();

            result.sf_squared =
                std::numeric_limits<double>::quiet_NaN();

        } else {

            result.sf =
                sum_sf[i] /
                static_cast<double>(count[i]);

            result.sf_squared =
                result.sf * result.sf;
        }

        results.push_back(result);
    }

    return SFResult(
        std::move(results)
    );
}

} // namespace esf