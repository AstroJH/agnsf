#pragma once

#include <vector>

#include <esf/sf_method.hpp>
#include <esf/lag_bins.hpp>
#include <esf/light_curve.hpp>
#include <esf/light_curve_view.hpp>
#include <esf/sf_result.hpp>

namespace esf {

/**
 * Ensemble SF computed by pooling all pairs from all light curves.
 *
 * For each lag bin, pair contributions from all light curves are
 * accumulated into the same BinAccumulator. Thus, the resulting
 * SF is weighted by the number of contributing pairs rather than
 * by the number of light curves.
 *
 * For a lag bin containing N pairs and the default estimator:
 *
 *   SF^2 =
 *       [sum(delta^2) - sum(error_i^2 + error_j^2)] / N
 *
 * The estimator (SFMethod) selects how SF^2 is formed from the
 * pooled pair statistics.
 *
 * No pair-level data are stored; only the accumulated statistics
 * required to compute the final SF are retained.
 */
class PooledESFCalculator {
public:
    SFResult calculate(
        const std::vector<LightCurve>& data,
        const LagBins& bins,
        SFMethod method = SFMethod::SecondOrder
    ) const;

    SFResult calculate(
        const std::vector<LightCurveView>& data,
        const LagBins& bins,
        SFMethod method = SFMethod::SecondOrder
    ) const;
};

} // namespace esf
