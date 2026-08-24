#pragma once

#include <esf/sf_result.hpp>
#include <esf/lag_bins.hpp>
#include <esf/light_curve.hpp>
#include <esf/light_curve_view.hpp>

namespace esf {

class SFCalculator {
public:
    SFResult calculate(
        const LightCurve& data,
        const LagBins& bins
    ) const;

    SFResult calculate(
        const LightCurveView& data,
        const LagBins& bins
    ) const;
};

} // namespace esf