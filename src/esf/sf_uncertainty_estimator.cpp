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
 * Map an uncertainty on SF^2 onto SF via the delta method (first order):
 *
 *   sigma_sf ~= sigma_sf2 / (2 * sf)
 *
 * Valid only when sf > 0; diverges as sf -> 0, so unestimated (NaN)
 * intervals are returned for non-positive sf.
 */
SFUncertainty delta_to_sf(
    double sf2,
    double sf,
    double sigma_sf2
)
{
    SFUncertainty result;

    if (!std::isfinite(sf2) || sf2 <= 0.0 ||
        !std::isfinite(sf) || sf <= 0.0 ||
        !std::isfinite(sigma_sf2)) {
        return result;
    }

    const double sigma_sf =
        sigma_sf2 / (2.0 * sf);

    result.lower = sf - sigma_sf;
    result.upper = sf + sigma_sf;

    return result;
}


/**
 * Map an uncertainty interval on SF^2 onto SF by propagating each
 * bound through sf = sqrt(SF^2) (lower clamped at 0). This gives an
 * asymmetric interval on sf and is used for the naive within-bin
 * statistical uncertainty.
 */
SFUncertainty sqrt_interval(
    double sf2,
    double delta_sf2
)
{
    SFUncertainty result;

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


// ---------------------------------------------------------------------
// Measurement uncertainty: propagation of the observation errors sigma_i.
// IMPLEMENTATION CONSTRAINT: pairs are treated as independent
// (shared-observation covariance is ignored).
// ---------------------------------------------------------------------

/**
 * Second-order family.
 *
 * Treat the underlying values V_i and measured values F_i as random
 * variables, related by
 *
 *   F_i = V_i + e_i,    e_i ~ N(0, sigma_i^2).
 *
 * For a pair p = (i, j), i != j, let
 *
 *   D_p = V_i - V_j
 *   d_p = F_i - F_j = D_p + (e_i - e_j)
 *
 * denote the underlying and measured differences, respectively. Define
 *
 *   s_p^2 = sigma_i^2 + sigma_j^2.
 *
 * Assuming independent measurement errors e_i and e_j,
 *
 *   Var(e_i-e_j) = s_p^2.
 * 
 * The measurement-induced variance of d_p^2 is
 *
 *   Var(d_p^2 | D_p)
 *        = Var(D_p^2 + 2 D_p(e_i-e_j) + (e_i-e_j)^2 | D_p)
 *        = 0 + 4 D_p^2 Var(e_i-e_j) + Var[(e_i-e_j)^2]
 *        = 4 D_p^2 s_p^2 + 2 s_p^4.
 *
 * This variance is *unchanged* for the bias-corrected form
 *
 *   X_p = d_p^2 - s_p^2,
 *
 * since the subtracted term is deterministic.
 *
 * Assuming independent pairs (WARNING: generally not exact),
 *
 *   Var_meas(SF^2) = 1/N^2 * sum_p[ Var(d_p^2 | D_p) ]
 *                  = 1/N^2 * sum_p[ 4 D_p^2 s_p^2 + 2 s_p^4 ].
 * 
 * In practice, D_p cannot be directly observed; d_p^2 is used to
 * approximate for D_p^2:
 *
 *   Var_meas(SF^2) = 1/N^2 * sum_p[ 4 d_p^2 s_p^2 + 2 s_p^4 ].
 */
class SecondOrderMeasurementUncertaintyEstimator
    : public SFMeasurementUncertaintyEstimator {
public:
    explicit SecondOrderMeasurementUncertaintyEstimator(
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

        if (n == 0) {
            return SFUncertainty{};
        }

        const double count = static_cast<double>(n);

        // sum_p[ 4 d_p^2 s_p^2 + 2 s_p^4 ]
        const double sum_measurement_var =
            4.0 * stats.sum_delta2_noise() +
            2.0 * stats.sum_noise2();

        const double sigma_sf2 =
            std::sqrt(sum_measurement_var) / count;

        const double sf2 = stats.sf_squared(method_);
        const double sf = stats.sf(method_);

        return delta_to_sf(sf2, sf, sigma_sf2);
    }

private:
    SFMethod method_;
};


 /**
  * Mean-absolute-deviation family.
  *
  * Treat the underlying values V_i and measured values F_i as random
  * variables, related by
  *
  *   F_i = V_i + e_i,    e_i ~ N(0, sigma_i^2).
  *
  * For a pair p = (i, j), i != j, let
  *
  *   D_p = V_i - V_j
  *   d_p = F_i - F_j = D_p + (e_i - e_j),
  *
  * and define
  *
  *   s_p^2 = sigma_i^2 + sigma_j^2.
  *
  * Assuming independent measurement errors e_i and e_j,
  *
  *   Var(e_i-e_j) = s_p^2.
  * 
  * SF^2 = pi/2 * <|d_p|>^2  (- <s_p^2>).
  * 
  * For small measurement errors relative to the underlying differences,
  * linear propagation through |d_p| gives
  *
  *   |d_p| = |D_p + (e_i - e_j)| ~= |D_p| + sign(D_p) * (e_i - e_j),
  *
  *   => Var(|d_p| | D_p) ~= s_p^2.
  * 
  * Assuming independent pairs (WARNING: generally not exact),
  * 
  *   Var_meas(<|d_p|>) ~= 1/N^2 * sum_p[ s_p^2 ].
  *
  * Propagating this uncertainty to SF^2 gives
  *
  *   Var_meas(SF^2) ~= pi^2 * <|d_p|>^2 * Var_meas(<|d_p|>).
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

        if (n == 0) {
            return SFUncertainty{};
        }

        const double count = static_cast<double>(n);

        const double mean_abs =
            stats.sum_abs_delta() / count;

        // pi^2 * <|d_p|>^2 * sum_p[ s_p^2 ] / N^2
        const double var_sf2 =
            kPi * kPi * // pi^2
            mean_abs * mean_abs * // <|d_p|>^2
            stats.sum_noise() / // sum_p[ s_p^2 ]
            (count * count); // N^2

        const double sigma_sf2 =
            std::sqrt(var_sf2);

        const double sf2 = stats.sf_squared(method_);
        const double sf = stats.sf(method_);

        return delta_to_sf(sf2, sf, sigma_sf2);
    }

private:
    SFMethod method_;
};


// ---------------------------------------------------------------------
// Naive within-bin statistical uncertainty: s_X / sqrt(N_pair).
// ---------------------------------------------------------------------

/**
 * Second-order family.
 *
 *   X_p = d_p^2 - s_p^2   (noise-corrected)
 *   X_p = d_p^2           (no-noise)
 *
 *   SE = sample_std(X_p) / sqrt(N).
 */
class SecondOrderWithinUncertaintyEstimator
    : public SFWithinUncertaintyEstimator {
public:
    explicit SecondOrderWithinUncertaintyEstimator(
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
 *   m = <|d_p|>
 *   SE(m) = sample_std(|d_p|) / sqrt(N)
 *   d(SF^2) = pi * m * SE(m).
 */
class MeanAbsoluteDeviationWithinUncertaintyEstimator
    : public SFWithinUncertaintyEstimator {
public:
    explicit MeanAbsoluteDeviationWithinUncertaintyEstimator(
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

        double sample_var =
            (
                stats.sum_delta_squared() -
                count * mean_abs * mean_abs
            ) / (count - 1.0);

        sample_var =
            std::max(sample_var, 0.0);

        const double se_abs =
            std::sqrt(sample_var / count);

        const double delta_sf2 =
            kPi * mean_abs * se_abs;

        const double sf2 =
            stats.sf_squared(method_);

        return sqrt_interval(sf2, delta_sf2);
    }

private:
    SFMethod method_;
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

    // Unknown SFMethod: no estimator -> unestimated (NaN).
    return nullptr;
}


std::unique_ptr<SFWithinUncertaintyEstimator>
make_sf_within_uncertainty_estimator(
    SFMethod method
)
{
    switch (method) {

        case SFMethod::SecondOrder:
        case SFMethod::SecondOrderNoNoise:
            return std::make_unique<
                SecondOrderWithinUncertaintyEstimator
            >(method);

        case SFMethod::MeanAbsoluteDeviation:
        case SFMethod::MeanAbsoluteDeviationNoNoise:
            return std::make_unique<
                MeanAbsoluteDeviationWithinUncertaintyEstimator
            >(method);
    }

    return nullptr;
}

} // namespace esf
} // namespace agnsf
