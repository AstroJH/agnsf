#pragma once

#include <vector>

#include <esf/sf_method.hpp>
#include <esf/sf_result.hpp>
#include <esf/sf_uncertainty.hpp>
#include <esf/lag_bins.hpp>
#include <core/light_curve.hpp>

namespace agnsf {
namespace esf {

/**
 * Ensemble SF obtained by first computing an SF for each light
 * curve and then combining the individual SF values.
 *
 * The per-curve SF is computed with the chosen estimator
 * (SFMethod). Two combination methods are then supported:
 *
 *   SqrtMeanSquared (default):
 *
 *       ESF(tau) = sqrt( <SF_k^2(tau)>_k )
 *
 *     For each lag bin, light curves with a finite SF^2
 *     contribute to the mean.
 *
 *   MeanSf:
 *
 *       ESF(tau) = <SF_k(tau)>_k
 *
 *     For each lag bin, light curves with a finite SF
 *     contribute to the mean.
 *
 * Each contributing light curve is weighted equally within
 * each lag bin.
 *
 * Uncertainty (UncertaintyConfig):
 *
 *   - measurement: per-curve measurement uncertainty (Analytic)
 *     propagated through the aggregation;
 *   - sampling: source-to-source scatter, either Analytic
 *     (std/sqrt(n) of the per-curve values) or curve-level
 *     Jackknife / Bootstrap.
 */
class SFEnsembleCalculator {
public:
    enum class Method {
        // ESF(tau) = sqrt( <SF_k^2(tau)>_k )
        SqrtMeanSquared,

        // ESF(tau) = <SF_k(tau)>_k
        MeanSf
    };

    SFResult calculate(
        const std::vector<agnsf::LightCurve>& data,
        const LagBins& bins,
        SFMethod sf_method = SFMethod::SecondOrder,
        Method method = Method::SqrtMeanSquared,
        const UncertaintyConfig& config = {},
        double redshift = 0.0
    ) const;

    SFResult calculate(
        const std::vector<agnsf::LightCurveView>& data,
        const LagBins& bins,
        SFMethod sf_method = SFMethod::SecondOrder,
        Method method = Method::SqrtMeanSquared,
        const UncertaintyConfig& config = {},
        double redshift = 0.0
    ) const;

    // Per-curve redshifts: one value per light curve.
    SFResult calculate(
        const std::vector<agnsf::LightCurve>& data,
        const LagBins& bins,
        SFMethod sf_method,
        Method method,
        const UncertaintyConfig& config,
        const std::vector<double>& redshifts
    ) const;

    SFResult calculate(
        const std::vector<agnsf::LightCurveView>& data,
        const LagBins& bins,
        SFMethod sf_method,
        Method method,
        const UncertaintyConfig& config,
        const std::vector<double>& redshifts
    ) const;
};

} // namespace esf
} // namespace agnsf
