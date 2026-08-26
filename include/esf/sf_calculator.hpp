#pragma once

#include <esf/sf_method.hpp>
#include <esf/sf_result.hpp>
#include <esf/sf_uncertainty.hpp>
#include <esf/lag_bins.hpp>
#include <core/light_curve.hpp>

namespace agnsf {
namespace esf {

/**
 * Structure function of a single light curve.
 *
 * The estimator (SFMethod) selects how SF^2 is formed from the
 * pair statistics in each lag bin. When the measurement uncertainty
 * is enabled (UncertaintyConfig::measurement == Analytic) each bin's
 * SFBinResult::measurement is filled from the within-bin standard
 * error of the mean propagated to sf.
 *
 * Sampling uncertainty is not defined for a single light curve and
 * is rejected by the config validation.
 */
class SFCalculator {
public:
    SFResult calculate(
        const agnsf::LightCurve& data,
        const LagBins& bins,
        SFMethod method = SFMethod::SecondOrder,
        const UncertaintyConfig& config = {}
    ) const;

    SFResult calculate(
        const agnsf::LightCurveView& data,
        const LagBins& bins,
        SFMethod method = SFMethod::SecondOrder,
        const UncertaintyConfig& config = {}
    ) const;
};

} // namespace esf
} // namespace agnsf
