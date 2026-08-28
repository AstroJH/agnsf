#pragma once

#include <memory>

#include <esf/bin_accumulator.hpp>
#include <esf/sf_method.hpp>
#include <core/uncertainty.hpp>

namespace agnsf {
namespace esf {

/**
 * Estimates the measurement uncertainty of a single-bin SF.
 *
 * "Measurement uncertainty" here means the uncertainty propagated from
 * the per-observation measurement errors sigma_i (e.g. photometric
 * errors) into the SF estimate; it depends explicitly on sigma_i.
 *
 * IMPLEMENTATION CONSTRAINT: the current closed-form propagation
 * treats the measurement contributions of different pairs as
 * independent (pairs that share an observation have zero covariance in
 * this approximation). A future, more rigorous model may account for
 * shared-observation covariance.
 *
 * Each SFMethod can provide its own propagation rule.
 */
class SFMeasurementUncertaintyEstimator {
public:
    virtual ~SFMeasurementUncertaintyEstimator() = default;

    /**
     * Estimate the measurement uncertainty of sf for one lag bin from
     * the bin's accumulated pair statistics.
     *
     * Returns an unestimated (NaN) interval when sf is not finite
     * (e.g. noise-dominated bins with negative SF^2).
     */
    virtual Uncertainty estimate(
        const BinAccumulator& stats
    ) const = 0;
};


/**
 * Estimates the naive within-bin statistical uncertainty of a
 * single-bin SF.
 *
 * This is the standard error of the per-pair mean under the
 * pair-independence approximation:
 *
 *   SE = sample_std(X) / sqrt(N_pair)
 *
 * where X is the per-pair statistic. Under a Gaussian process this is
 * a lower-bound-type approximation. It is NOT measurement error.
 */
class SFWithinUncertaintyEstimator {
public:
    virtual ~SFWithinUncertaintyEstimator() = default;

    virtual Uncertainty estimate(
        const BinAccumulator& stats
    ) const = 0;
};


/**
 * Create the measurement-uncertainty estimator (sigma_i propagation)
 * for an SFMethod. Returns nullptr for SFMethods without an
 * estimator, which callers treat as unestimated (NaN).
 */
std::unique_ptr<SFMeasurementUncertaintyEstimator>
make_sf_measurement_uncertainty_estimator(
    SFMethod method
);


/**
 * Create the naive within-bin statistical-uncertainty estimator for
 * an SFMethod. Returns nullptr for SFMethods without an estimator.
 */
std::unique_ptr<SFWithinUncertaintyEstimator>
make_sf_within_uncertainty_estimator(
    SFMethod method
);

} // namespace esf
} // namespace agnsf
