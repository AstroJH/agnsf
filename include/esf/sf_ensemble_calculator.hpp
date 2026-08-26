#pragma once

#include <vector>

#include <esf/sf_method.hpp>
#include <esf/lag_bins.hpp>
#include <esf/light_curve.hpp>
#include <esf/light_curve_view.hpp>
#include <esf/sf_result.hpp>

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
        const std::vector<LightCurve>& data,
        const LagBins& bins,
        SFMethod sf_method = SFMethod::SecondOrder,
        Method method = Method::SqrtMeanSquared
    ) const;

    SFResult calculate(
        const std::vector<LightCurveView>& data,
        const LagBins& bins,
        SFMethod sf_method = SFMethod::SecondOrder,
        Method method = Method::SqrtMeanSquared
    ) const;
};

} // namespace esf
