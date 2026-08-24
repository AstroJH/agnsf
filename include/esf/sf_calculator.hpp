#pragma once

#include <esf/sf_result.hpp>
#include <esf/lag_bins.hpp>
#include <esf/light_curve.hpp>

namespace esf {

class SFCalculator {
public:
    SFResult calculate(
        const LightCurve& data,
        const LagBins& bins
    ) const;
};

} // namespace esf