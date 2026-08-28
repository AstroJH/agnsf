#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>

namespace agnsf {

/**
 * Uncertainty interval attached to the structure function value `sf`.
 *
 * The interval is stored as absolute bounds on sf:
 *
 *   lower <= sf <= upper
 *
 * NaN means "not estimated". A symmetric estimate is represented by
 * lower == upper (i.e. sf +/- sigma). Asymmetric intervals are
 * supported directly.
 */
struct Uncertainty {
    double lower =
        std::numeric_limits<double>::quiet_NaN();

    double upper =
        std::numeric_limits<double>::quiet_NaN();

    /**
     * True when at least one side of the interval has been estimated.
     */
    bool estimated() const noexcept
    {
        const bool lower_finite = (lower == lower); // NaN check
        const bool upper_finite = (upper == upper);

        return lower_finite || upper_finite;
    }
};


/**
 * Uncertainty estimation method.
 *
 * `Off` disables estimation. The other values select how a given
 * uncertainty component is computed. The set of valid methods depends
 * on the component and the calculator:
 *
 *   measurement (SF or ESF) — uncertainty propagated from the
 *     per-observation measurement errors sigma_i (e.g. photometric
 *     errors):
 *     - Analytic: closed-form error propagation.
 *       IMPLEMENTATION CONSTRAINT: the measurement contributions of different pairs
 *       are treated as independent; pairs that share an observation
 *       therefore have zero covariance in this approximation. A
 *       future, more rigorous model may account for shared-observation
 *       covariance.
 *     - MonteCarlo: observation-level perturbation, f_i* = f_i + e_i
 *       with e_i ~ N(0, sigma_i^2), which naturally preserves the
 *       covariance induced by shared observations.
 *
 *   within (SF or ESF) — naive statistical uncertainty of the bin
 *     mean: sample_std(X) / sqrt(N_pair) over the per-pair statistics
 *     X, treating pairs as independent. This is a lower-bound-type
 *     approximation (e.g. under a Gaussian process) and is NOT
 *     measurement error.
 *
 *   sampling (ESF only, source-to-source):
 *     - Analytic: std of per-curve values / sqrt(n) (aggregated ESF).
 *     - Jackknife / Bootstrap: curve-level resampling.
 */
enum class UncertaintyMethod {
    Off = 0,
    Analytic,
    MonteCarlo,
    Jackknife,
    Bootstrap
};


/**
 * User configuration for uncertainty estimation.
 *
 * The default (all Off) reproduces the original point-estimate
 * behaviour exactly.
 *
 * `n_bootstrap` is the number of resamples for Bootstrap sampling and
 * the number of realizations for MonteCarlo measurement.
 * `bootstrap_seed` fixes the RNG so results are reproducible.
 */
struct UncertaintyConfig {
    // Measurement uncertainty propagated from sigma_i.
    UncertaintyMethod measurement = UncertaintyMethod::Off;

    // Naive within-bin statistical uncertainty (s_X / sqrt(N_pair)).
    UncertaintyMethod within = UncertaintyMethod::Off;

    // Source-to-source sampling uncertainty (ESF only).
    UncertaintyMethod sampling = UncertaintyMethod::Off;

    // Number of bootstrap resamples / Monte Carlo realizations.
    std::size_t n_bootstrap = 100;

    // RNG seed for resampling / perturbation; 0 gives a fixed default.
    std::uint32_t bootstrap_seed = 0;
};

} // namespace agnsf
