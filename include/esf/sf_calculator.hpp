#pragma once

#include <esf/sf_method.hpp>
#include <esf/sf_result.hpp>
#include <esf/lag_bins.hpp>
#include <core/light_curve.hpp>

namespace agnsf {
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
        const agnsf::LightCurve& data,
        const LagBins& bins,
        SFMethod method = SFMethod::SecondOrder
    ) const;

    SFResult calculate(
        const agnsf::LightCurveView& data,
        const LagBins& bins,
        SFMethod method = SFMethod::SecondOrder
    ) const;
};

} // namespace esf
} // namespace agnsf
