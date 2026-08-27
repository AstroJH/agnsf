#include <esf/sf_uncertainty_estimator.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

namespace agnsf {
namespace esf {

namespace {

constexpr double kPi = 3.14159265358979323846;


/**
 * Map an uncertainty on SF^2 (working scale) onto sf.
 *
 * The interval [sf2 - d, sf2 + d] is propagated through
 * sf = sqrt(sf2); the lower bound is clamped at 0 because sf is
 * non-negative by construction. The result is asymmetric whenever the
 * interval on SF^2 is symmetric, which is why asymmetric intervals are
 * supported from the start.
 */
SFUncertainty sqrt_interval(
    double sf2,
    double delta_sf2
)
{
    SFUncertainty result;

    // sf = sqrt(sf^2) is NaN for negative SF^2 (noise-dominated bin),
    // so no meaningful interval can be attached to sf.
    if (!std::isfinite(sf2) || sf2 < 0.0 ||
        !std::isfinite(delta_sf2)) {
        return result;
    }

    result.lower =
        std::sqrt(
            std::max(sf2 - delta_sf2, 0.0)
        );

    result.upper =
        std::sqrt(sf2 + delta_sf2);

    return result;
}


/**
 * Second-order family.
 *
 * SF^2 is the mean of the per-pair values
 *
 *   x_k = delta_k^2 - (sigma_i^2 + sigma_j^2)     (noise-corrected)
 *   x_k = delta_k^2                               (no-noise)
 *
 * The measurement uncertainty is the standard error of that mean:
 *
 *   se = sample_std(x) / sqrt(N)
 *
 * where sample_std uses the accumulated moments. The interval on
 * SF^2 is then propagated through sf = sqrt(SF^2).
 *
 * NOTE: pairs sharing points are treated as independent; this is an
 * approximation chosen for the initial implementation (the future
 * Monte Carlo perturbation option will capture pair correlations).
 */
class SecondOrderMeasurementUncertaintyEstimator
    : public SFMeasurementUncertaintyEstimator {
public:
    explicit SecondOrderMeasurementUncertaintyEstimator(
        SFMethod method
    )
        : method_(method),
          subtract_noise_(method == SFMethod::SecondOrder)
    {
    }

    SFUncertainty estimate(
        const BinAccumulator& stats
    ) const override
    {
        const std::size_t n = stats.count();

        if (n < 2) {
            return SFUncertainty{};
        }

        const double count = static_cast<double>(n);

        // Mean and second moment of the per-pair values x_k.
        const double sum_x =
            subtract_noise_
                ? stats.sum_delta_squared() -
                  stats.sum_noise()
                : stats.sum_delta_squared();

        const double sum_x2 =
            subtract_noise_
                ? stats.sum_delta4() -
                  2.0 * stats.sum_delta2_noise() +
                  stats.sum_noise2()
                : stats.sum_delta4();

        const double mean = sum_x / count;

        // Sample variance of x_k; clamp tiny negative rounding.
        double sample_var =
            (sum_x2 - count * mean * mean) /
            (count - 1.0);

        sample_var =
            std::max(sample_var, 0.0);

        const double se =
            std::sqrt(sample_var / count);

        const double sf2 =
            stats.sf_squared(method_);

        return sqrt_interval(sf2, se);
    }

private:
    SFMethod method_;
    bool subtract_noise_;
};


/**
 * Mean-absolute-deviation family.
 *
 * SF^2 = pi/2 * <|delta|>^2  (- <noise>).
 *
 * The measurement uncertainty is obtained from the standard error of
 * the mean absolute difference and propagated through the non-linear
 * factor pi/2 * x^2 (delta method):
 *
 *   se_abs   = sample_std(|delta|) / sqrt(N)
 *   d(SF^2)  = pi * <|delta|> * se_abs
 *
 * The measurement-noise term is treated as known exactly.
 */
class MeanAbsoluteDeviationMeasurementUncertaintyEstimator
    : public SFMeasurementUncertaintyEstimator {
public:
    explicit MeanAbsoluteDeviationMeasurementUncertaintyEstimator(
        SFMethod method
    )
        : method_(method)
    {
    }

    SFUncertainty estimate(
        const BinAccumulator& stats
    ) const override
    {
        const std::size_t n = stats.count();

        if (n < 2) {
            return SFUncertainty{};
        }

        const double count = static_cast<double>(n);

        const double mean_abs =
            stats.sum_abs_delta() / count;

        // Sample variance of |delta| (sum(|delta|^2) == sum(delta^2));
        // clamp tiny negative rounding.
        double sample_var =
            (
                stats.sum_delta_squared() -
                count * mean_abs * mean_abs
            ) / (count - 1.0);

        sample_var =
            std::max(sample_var, 0.0);

        const double se_abs =
            std::sqrt(sample_var / count);

        // d/dx (pi/2 * x^2) = pi * x.
        const double delta_sf2 =
            kPi * mean_abs * se_abs;

        const double sf2 =
            stats.sf_squared(method_);

        return sqrt_interval(sf2, delta_sf2);
    }

private:
    SFMethod method_;
};


/**
 * Fallback used for SFMethods without a dedicated estimator.
 *
 * Returns an unestimated interval so that adding a new SFMethod never
 * requires every caller to handle a missing estimator.
 */
class NoOpMeasurementUncertaintyEstimator
    : public SFMeasurementUncertaintyEstimator {
public:
    SFUncertainty estimate(
        const BinAccumulator&
    ) const override
    {
        return SFUncertainty{};
    }
};

} // namespace


std::unique_ptr<SFMeasurementUncertaintyEstimator>
make_sf_measurement_uncertainty_estimator(
    SFMethod method
)
{
    switch (method) {

        case SFMethod::SecondOrder:
        case SFMethod::SecondOrderNoNoise:
            return std::make_unique<
                SecondOrderMeasurementUncertaintyEstimator
            >(method);

        case SFMethod::MeanAbsoluteDeviation:
        case SFMethod::MeanAbsoluteDeviationNoNoise:
            return std::make_unique<
                MeanAbsoluteDeviationMeasurementUncertaintyEstimator
            >(method);
    }

    return std::make_unique<
        NoOpMeasurementUncertaintyEstimator
    >();
}

} // namespace esf
} // namespace agnsf
