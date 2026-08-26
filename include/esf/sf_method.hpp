#pragma once

namespace agnsf {
namespace esf {

/**
 * Structure-function estimator.
 *
 * Every estimator is defined per lag bin from the pair statistics
 * accumulated by BinAccumulator:
 *
 *   SecondOrder (default):
 *       SF^2 = <delta^2> - <sigma_i^2 + sigma_j^2>
 *
 *   SecondOrderNoNoise:
 *       SF^2 = <delta^2>
 *
 *   MeanAbsoluteDeviation:
 *       SF^2 = pi/2 * <|delta|>^2 - <sigma_i^2 + sigma_j^2>
 *
 *   MeanAbsoluteDeviationNoNoise:
 *       SF^2 = pi/2 * <|delta|>^2
 *
 * where delta = value_j - value_i and the averages run over all
 * pairs in the lag bin. In every case SF = sqrt(SF^2) when SF^2 is
 * finite and non-negative, otherwise SF is NaN.
 */
enum class SFMethod {
    SecondOrder,
    SecondOrderNoNoise,
    MeanAbsoluteDeviation,
    MeanAbsoluteDeviationNoNoise
};

} // namespace esf
} // namespace agnsf
