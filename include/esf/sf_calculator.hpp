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
 * The estimator (SFMethod) determines how SF^2 is calculated from
 * pair statistics in each lag bin.
 *
 * Uncertainty estimation is controlled by UncertaintyConfig.
 */
class SFCalculator {
public:

    /**
     * Calculate the structure function of a single light curve.
     *
     * @param data      Input light curve.
     * @param bins      Lag bins.
     * @param method    Structure-function estimator.
     * @param config    Uncertainty configuration.
     * @param redshift  Source redshift used for rest-frame lag correction.
     */
    SFResult calculate(
        const agnsf::LightCurve& data,
        const LagBins& bins,
        SFMethod method = SFMethod::SecondOrder,
        const UncertaintyConfig& config = {},
        double redshift = 0.0
    ) const;

    SFResult calculate(
        const agnsf::LightCurveView& data,
        const LagBins& bins,
        SFMethod method = SFMethod::SecondOrder,
        const UncertaintyConfig& config = {},
        double redshift = 0.0
    ) const;
};

} // namespace esf
} // namespace agnsf
