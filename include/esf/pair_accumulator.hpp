#pragma once

#include <vector>

#include <esf/bin_accumulator.hpp>
#include <esf/lag_bins.hpp>
#include <esf/light_curve.hpp>

namespace esf {

/**
 * Process one light curve and accumulate all of its pair
 * contributions into one local accumulator array.
 *
 * Input:
 *
 *   light_curve
 *   bins
 *
 * Output:
 *
 *   one BinAccumulator for each lag bin.
 *
 * No SF values are calculated here. The returned accumulators
 * contain only the pair-level statistics required for the
 * later Pooled ESF calculation.
 */
std::vector<BinAccumulator> accumulate_light_curve(
    const LightCurve& light_curve,
    const LagBins& bins
);

}