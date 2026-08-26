#pragma once

#include <esf/sf_method.hpp>
#include <esf/sf_result.hpp>
#include <esf/lag_bins.hpp>
#include <esf/light_curve.hpp>
#include <esf/light_curve_view.hpp>

namespace esf {

/**
 * Structure function of a single light curve.
 *
 * The estimator (SFMethod) selects how SF^2 is formed from the
 * pair statistics in each lag bin.
 */
class SFCalculator {
public:
    SFResult calculate(
        const LightCurve& data,
        const LagBins& bins,
        SFMethod method = SFMethod::SecondOrder
    ) const;

    SFResult calculate(
        const LightCurveView& data,
        const LagBins& bins,
        SFMethod method = SFMethod::SecondOrder
    ) const;
};

} // namespace esf
