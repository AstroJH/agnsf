#pragma once

#include <cstddef>

namespace esf {

/**
 * Accumulates pair-level statistics for a single lag bin.
 *
 * Each pair contributes:
 *
 *   delta^2
 *
 * and the corresponding measurement-noise term:
 *
 *   error_i^2 + error_j^2.
 *
 * The raw pair data are not stored. Instead, only the accumulated
 * statistics required to compute the noise-corrected SF are retained.
 *
 * For N contributing pairs:
 *
 *   SF^2 =
 *       [sum(delta^2) - sum(error_i^2 + error_j^2)] / N
 *
 * The accumulator is intended to be used during a single SF
 * calculation and represents the statistics of one lag bin.
 */
class BinAccumulator {
public:
    BinAccumulator() noexcept = default;

    /**
     * Add one pair contribution to the accumulator.
     *
     * @param delta   Difference in the observed quantity between
     *                the two measurements.
     * @param error_i Measurement uncertainty of the first point.
     * @param error_j Measurement uncertainty of the second point.
     */
    void add(
        double delta,
        double error_i,
        double error_j
    ) noexcept;

    /**
     * Merge the accumulated statistics from another accumulator.
     *
     * The two accumulators must correspond to the same lag bin.
     */
    void merge(
        const BinAccumulator& other
    ) noexcept;

    /**
     * Number of pairs accumulated in this bin.
     */
    std::size_t count() const noexcept;

    /**
     * Sum of squared pair differences:
     *
     *   sum(delta^2)
     */
    double sum_delta_squared() const noexcept;

    /**
     * Sum of the pair-wise measurement-noise terms:
     *
     *   sum(error_i^2 + error_j^2)
     */
    double sum_noise() const noexcept;

    /**
     * Noise-corrected second-order structure function squared:
     *
     *   SF^2 =
     *       [sum(delta^2) - sum_noise] / count
     *
     * Returns NaN if the bin contains no pairs.
     */
    double sf_squared() const noexcept;

    /**
     * Structure function:
     *
     *   SF = sqrt(SF^2)
     *
     * Returns NaN if SF^2 is not finite or is negative.
     */
    double sf() const noexcept;

private:
    // Number of pair contributions.
    std::size_t count_ = 0;

    // Accumulated squared pair differences.
    double sum_delta_squared_ = 0.0;

    // Accumulated measurement-noise contribution.
    double sum_noise_ = 0.0;
};

} // namespace esf
