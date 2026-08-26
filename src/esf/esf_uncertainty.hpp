// Internal helpers shared by the ESF calculators.

#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <esf/sf_uncertainty.hpp>

namespace agnsf {
namespace esf {
namespace detail {

/**
 * Analytic sampling interval on the given per-curve values:
 *
 *   mean +/- sample_std(values) / sqrt(n)
 *
 * Returns an unestimated interval when n < 2.
 */
SFUncertainty analytic_interval(
    const std::vector<double>& values
);


/**
 * Curve-level jackknife interval.
 *
 * For each k the statistic is recomputed from all curves except k and
 * the jackknife variance is
 *
 *   var = (n - 1) / n * sum_k ( theta_{(-k)} - mean(theta_{(-k)}) )^2
 *
 * The interval is centered on the full-sample statistic. `sqrt_mean`
 * selects whether the per-replicate statistic is the mean itself
 * (MeanSf) or sqrt(mean) (SqrtMeanSquared); negative arguments are
 * clamped to 0 so that sf stays real.
 *
 * Returns an unestimated interval when n < 2.
 */
SFUncertainty jackknife_interval(
    const std::vector<double>& values,
    bool sqrt_mean
);


/**
 * Curve-level bootstrap interval (percentile method, 16% / 84%).
 *
 * Resamples the per-curve values with replacement `n_bootstrap` times
 * and uses the 16th / 84th percentiles of the resampled statistic
 * (mean, or sqrt(mean) when `sqrt_mean`) as the interval. The 16/84
 * percentiles give a ~1-sigma central interval and are naturally
 * asymmetric. The same seed reproduces the same result.
 *
 * Returns an unestimated interval when n < 2 or when fewer than two
 * resampled statistics are finite.
 */
SFUncertainty bootstrap_interval(
    const std::vector<double>& values,
    std::size_t n_bootstrap,
    std::uint32_t seed,
    bool sqrt_mean
);


/**
 * Map an interval from the working scale (sf or sf^2) onto sf.
 *
 * When `sqrt_transform` is true each bound is propagated through
 * sf = sqrt(working) with the lower bound clamped at 0; otherwise the
 * interval is returned unchanged.
 */
SFUncertainty map_interval_to_sf(
    const SFUncertainty& interval,
    bool sqrt_transform
);

} // namespace detail
} // namespace esf
} // namespace agnsf
