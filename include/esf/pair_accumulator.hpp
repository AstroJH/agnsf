#pragma once

#include <vector>

#include <esf/bin_accumulator.hpp>
#include <esf/lag_bins.hpp>
#include <core/light_curve.hpp>

namespace agnsf {
namespace esf {

/**
 * Process one light curve and accumulate all of its pair
 * contributions into one local accumulator array.
 *
 * Input:
 *
 *   light_curve
 *   bins
 *   redshift (optional)
 *
 * Output:
 *
 *   one BinAccumulator for each lag bin.
 *
 * When `redshift` is non-zero, lags are converted to the source rest
 * frame before binning:
 *
 *   dt_rest = dt_obs / (1 + z)
 *
 * `redshift` must be > -1.
 *
 * No SF values are calculated here. The returned accumulators
 * contain only the pair-level statistics required for the
 * later Pooled ESF calculation.
 */
std::vector<BinAccumulator> accumulate_light_curve(
    const agnsf::LightCurve& light_curve,
    const LagBins& bins,
    double redshift = 0.0
);

std::vector<BinAccumulator> accumulate_light_curve(
    const agnsf::LightCurveView& light_curve,
    const LagBins& bins,
    double redshift = 0.0
);
} // namespace esf
} // namespace agnsf
