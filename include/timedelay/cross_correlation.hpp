#pragma once

#include <cstddef>
#include <limits>
#include <vector>

#include <core/light_curve.hpp>

namespace agnsf {
namespace timedelay {

/**
 * Cross-correlation method for time-delay estimation.
 *
 *   Dcf:  discrete correlation function (Edelson & Krolik 1988).
 *         Uses every pair of observations directly, binned by lag;
 *         no interpolation.
 *
 *   Iccf: interpolated cross-correlation function
 *         (Gaskell & Peterson 1987). The response light curve is
 *         linearly interpolated onto the continuum times shifted by
 *         each trial lag, and the Pearson correlation is evaluated on
 *         the overlapping points.
 */
enum class CrossCorrelationMethod {
    Dcf = 0,
    Iccf
};


/**
 * How to turn the cross-correlation curve into a single lag estimate.
 *
 *   Peak:     the lag maximizing the correlation.
 *   Centroid: the flux-weighted centroid of the correlation curve
 *             above `centroid_threshold * peak`.
 */
enum class LagEstimate {
    Peak = 0,
    Centroid
};


/**
 * Uniform grid of trial lags (in the same time units as the light
 * curves). Positive lags mean the response lags the continuum.
 */
struct LagGrid {
    double min = -50.0;
    double max = 50.0;
    double step = 1.0;

    /** Number of trial lags (0 if the grid is invalid). */
    std::size_t size() const noexcept;
};


/**
 * Options for cross_correlate().
 */
struct CrossCorrelationConfig {
    CrossCorrelationMethod method = CrossCorrelationMethod::Dcf;

    // DCF bin width in time units.
    double dcf_bin_width = 1.0;

    // Fraction of the peak used as the centroid threshold (0..1).
    double centroid_threshold = 0.8;

    // Minimum number of overlapping points for a valid correlation.
    std::size_t min_overlap = 3;
};


/**
 * Result of a cross-correlation lag search.
 *
 * `tau`, `ccf`, and `count` are the trial-lag grid, the correlation
 * values (NaN where the overlap is insufficient), and the number of
 * contributing points per bin (DCF: pairs; ICCF: overlapping points).
 */
struct LagResult {
    double lag_peak =
        std::numeric_limits<double>::quiet_NaN();

    double lag_centroid =
        std::numeric_limits<double>::quiet_NaN();

    double peak_value =
        std::numeric_limits<double>::quiet_NaN();

    std::vector<double> tau;
    std::vector<double> ccf;
    std::vector<std::size_t> count;
};


/**
 * Cross-correlate a continuum light curve against a response light
 * curve over a grid of trial lags.
 *
 * Convention: positive lag means the response is delayed relative to
 * the continuum (response(t + lag) correlates with continuum(t)).
 *
 * @throws std::invalid_argument if the lag grid is invalid or the
 *         light curves are empty.
 */
LagResult cross_correlate(
    const agnsf::LightCurve& continuum,
    const agnsf::LightCurve& response,
    const LagGrid& grid,
    const CrossCorrelationConfig& config = {}
);


/**
 * Extract the requested lag estimate (peak or centroid) from a result.
 */
double lag_value(
    const LagResult& result,
    LagEstimate estimate
);

} // namespace timedelay
} // namespace agnsf
