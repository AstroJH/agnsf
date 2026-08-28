#pragma once

#include <vector>

#include <esf/sf_method.hpp>
#include <esf/sf_result.hpp>
#include <core/uncertainty.hpp>
#include <esf/lag_bins.hpp>
#include <core/light_curve.hpp>

namespace agnsf {
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
 *
 * Uncertainty (UncertaintyConfig):
 *
 *   - measurement: analytic within-bin standard error on the pooled
 *     pair statistics;
 *   - sampling: source-to-source scatter via curve-level Jackknife or
 *     Bootstrap (Analytic is not defined for pooled ESF because no
 *     per-curve statistics are averaged).
 */
class PooledESFCalculator {
public:
    SFResult calculate(
        const std::vector<agnsf::LightCurve>& data,
        const LagBins& bins,
        SFMethod method = SFMethod::SecondOrder,
        const UncertaintyConfig& config = {},
        double redshift = 0.0
    ) const;

    SFResult calculate(
        const std::vector<agnsf::LightCurveView>& data,
        const LagBins& bins,
        SFMethod method = SFMethod::SecondOrder,
        const UncertaintyConfig& config = {},
        double redshift = 0.0
    ) const;

    // Per-curve redshifts: one value per light curve.
    SFResult calculate(
        const std::vector<agnsf::LightCurve>& data,
        const LagBins& bins,
        SFMethod method,
        const UncertaintyConfig& config,
        const std::vector<double>& redshifts
    ) const;

    SFResult calculate(
        const std::vector<agnsf::LightCurveView>& data,
        const LagBins& bins,
        SFMethod method,
        const UncertaintyConfig& config,
        const std::vector<double>& redshifts
    ) const;
};

} // namespace esf
} // namespace agnsf
