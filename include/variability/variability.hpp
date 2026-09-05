#pragma once

#include <cstddef>

#include <core/light_curve.hpp>
#include <core/uncertainty.hpp>

namespace agnsf {
namespace variability {

/**
 * Single-curve variability statistics.
 *
 * All measures are computed from one light curve (values + per-point
 * measurement errors).
 *
 * Notation:
 *
 *   n        number of points
 *   x_i      values, xbar = mean
 *   sigma_i  per-point errors
 *   eps2     mean square error = <sigma_i^2> (+ err_sys^2 if given)
 *   S2       sample variance = sum (x_i - xbar)^2 / (n - 1)
 *
 * Definitions follow common AGN variability practice
 * (Vaughan et al. 2003, MNRAS 345, 1271, and references therein):
 *
 *   sigma_m (intrinsic amplitude)     = sqrt(max(S2 - eps2, 0))
 *   sigma^2_NXS (norm. excess var.)   = (S2 - eps2) / xbar^2
 *   sigma^2_XS  (excess variance)     = S2 - eps2
 *   F_var       (frac. var. amplitude)= sqrt(max(sigma^2_NXS, 0))
 *
 * with the Vaughan et al. (2003) analytic uncertainties
 *
 *   err(F_var) = sqrt{ eps2/(n xbar^2)
 *                      + (1/(2n)) [eps2/(xbar^2 F_var)]^2 }
 *   err(sigma^2_NXS) = 2 F_var * err(F_var)
 *   err(sigma^2_XS)  = xbar^2 * err(sigma^2_NXS)
 *
 * Uncertainties are symmetric and stored as absolute bounds
 * [value - err, value + err] in an agnsf::Uncertainty.
 */
struct Options {
    // Systematic (floor) error added in quadrature to the per-point
    // errors when subtracting the noise contribution.
    double err_sys = 0.0;

    // Use the inverse-variance weighted mean (instead of the plain
    // mean) inside the variance / sigma_m estimates.
    bool weighted = false;
};


/**
 * All statistics of one light curve.
 *
 * `valid` is true when the data are sufficient (n >= 2, finite
 * values). Individual fields are NaN when not defined (e.g. zero
 * mean for F_var, n < 3 for the von Neumann ratio, no positive
 * errors for the weighted mean).
 */
struct Statistics {
    std::size_t n = 0;

    double mean = 0.0;

    double weighted_mean = 0.0;
    double weighted_mean_error = 0.0;

    double stddev = 0.0;
    double stddev_error = 0.0;

    double peak_to_peak = 0.0;
    double peak_to_peak_noise_corrected = 0.0;

    double sigma_m = 0.0;

    double fvar = 0.0;
    Uncertainty fvar_uncertainty;

    double nxs = 0.0;
    Uncertainty nxs_uncertainty;

    double xs = 0.0;
    Uncertainty xs_uncertainty;

    double chi2 = 0.0;
    double chi2_dof = 0.0;
    double chi2_q = 0.0;

    double von_neumann = 0.0;

    bool valid = false;
};


/**
 * Compute all variability statistics of one light curve.
 */
Statistics measure(
    const LightCurveView& curve,
    const Options& options = {}
);

} // namespace variability
} // namespace agnsf
