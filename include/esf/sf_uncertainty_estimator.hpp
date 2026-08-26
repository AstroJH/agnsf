#pragma once

#include <memory>

#include <esf/bin_accumulator.hpp>
#include <esf/sf_method.hpp>
#include <esf/sf_uncertainty.hpp>

namespace agnsf {
namespace esf {

/**
 * Interface for estimating the measurement uncertainty of a single-bin SF.
 *
 * Each SFMethod (or family of closely related methods) can provide its own
 * uncertainty estimator. The propagation formula is implemented by the
 * estimator rather than by the caller, allowing different SF methods to use
 * different uncertainty models.
 * 
 * Analytic propagation is the current implementation;
 * TODO: Monte Carlo / bootstrap based estimators can be
 * added behind the same interface later.
 */
class SFMeasurementUncertaintyEstimator {
public:
    virtual ~SFMeasurementUncertaintyEstimator() = default;

    /**
     * Estimate the measurement uncertainty of SF for one lag bin.
     *
     * @param stats Accumulated pair-level statistics for the bin.
     *
     * @return An uncertainty interval for SF. Returns an unestimated
     *         (NaN) interval if the bin contains fewer than two pairs
     *         or if the corresponding SF is not finite.
     */
    virtual SFUncertainty estimate(
        const BinAccumulator& stats
    ) const = 0;
};


/**
 * Factory function that creates the measurement-uncertainty estimator
 * corresponding to the specified SFMethod.
 *
 * The returned estimator encapsulates the uncertainty-propagation rule
 * for that SF method. Closely related SF methods may share the same
 * estimator implementation.
 *
 * If no estimator is currently implemented for a given SFMethod, the
 * factory returns an estimator that reports the uncertainty as
 * unestimated (NaN), allowing new SF methods to be added without
 * immediately requiring an uncertainty implementation.
 */
std::unique_ptr<SFMeasurementUncertaintyEstimator>
make_sf_measurement_uncertainty_estimator(
    SFMethod method
);

} // namespace esf
} // namespace agnsf
