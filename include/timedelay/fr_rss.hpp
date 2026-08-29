#pragma once

#include <cstddef>
#include <cstdint>

#include <core/light_curve.hpp>
#include <core/uncertainty.hpp>

#include <timedelay/cross_correlation.hpp>

namespace agnsf {
namespace timedelay {

/**
 * Options for the FR/RSS Monte Carlo lag uncertainty
 * (Peterson et al. 1998).
 *
 *   Flux randomization (FR): perturb every measurement by its error,
 *       f_i* = f_i + e_i, e_i ~ N(0, sigma_i^2).
 *
 *   Random subset selection (RSS): resample a random subset of the
 *       observed epochs (bootstrap with replacement) for each
 *       realization.
 *
 * Both are enabled by default. The 16th/84th percentiles of the
 * realization lag distribution form the reported interval.
 */
struct FRRSSConfig {
    std::size_t n_realizations = 1000;
    std::uint32_t seed = 0;
    bool flux_randomization = true;
    bool random_subset = true;
};


/**
 * FR/RSS uncertainty of the chosen lag estimate.
 *
 * Returns an unestimated (NaN) interval when fewer than two
 * realizations produce a finite lag.
 */
agnsf::Uncertainty lag_uncertainty(
    const agnsf::LightCurve& continuum,
    const agnsf::LightCurve& response,
    const LagGrid& grid,
    LagEstimate estimate,
    const CrossCorrelationConfig& ccf_config = {},
    const FRRSSConfig& config = {}
);

} // namespace timedelay
} // namespace agnsf
