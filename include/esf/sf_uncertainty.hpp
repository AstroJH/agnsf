#pragma once

#include <cstddef>
#include <cstdint>
#include <limits>

namespace agnsf {
namespace esf {

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
struct SFUncertainty {
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
 *   measurement (SF or ESF):
 *     - Analytic (default when enabled): within-bin standard error of
 *       the mean, propagated to SF. TODO: Monte Carlo perturbation may
 *       be added here in the future.
 *
 *   sampling (ESF only, source-to-source):
 *     - Analytic: std of per-curve values / sqrt(n) (aggregated ESF).
 *     - Jackknife / Bootstrap: curve-level resampling.
 */
enum class UncertaintyMethod {
    Off = 0,
    Analytic,
    Jackknife,
    Bootstrap
};


/**
 * User configuration for uncertainty estimation.
 *
 * The default (all Off) reproduces the current behaviour exactly.
 *
 * `n_bootstrap` and `bootstrap_seed` are used only when `sampling` is
 * Bootstrap. A fixed seed makes the result reproducible.
 */
struct UncertaintyConfig {
    // Measurement uncertainty (per-bin analytic propagation).
    UncertaintyMethod measurement = UncertaintyMethod::Off;

    // Source-to-source sampling uncertainty (ESF only).
    UncertaintyMethod sampling = UncertaintyMethod::Off;

    // Number of bootstrap resamples (sampling == Bootstrap).
    std::size_t n_bootstrap = 100;

    // RNG seed for bootstrap resampling; 0 gives a fixed default.
    std::uint32_t bootstrap_seed = 0;
};

} // namespace esf
} // namespace agnsf
