#include <cmath>
#include <limits>

#include <esf/bin_accumulator.hpp>

namespace agnsf {
namespace esf {

namespace {

constexpr double kPi = 3.14159265358979323846;

} // namespace


void BinAccumulator::add(
    double delta,
    double error_i,
    double error_j
) noexcept
{
    const double delta_squared = delta * delta;
    const double noise =
        error_i * error_i +
        error_j * error_j;

    sum_delta_squared_ += delta_squared;
    sum_abs_delta_ += std::abs(delta);
    sum_noise_ += noise;

    // Additional pair-level sums used by the analytic
    // measurement-uncertainty estimators.
    sum_delta4_ += delta_squared * delta_squared;
    sum_delta2_noise_ += delta_squared * noise;
    sum_noise2_ += noise * noise;

    ++count_;
}


void BinAccumulator::merge(
    const BinAccumulator& other
) noexcept
{
    count_ += other.count_;
    sum_delta_squared_ += other.sum_delta_squared_;
    sum_abs_delta_ += other.sum_abs_delta_;
    sum_noise_ += other.sum_noise_;
    sum_delta4_ += other.sum_delta4_;
    sum_delta2_noise_ += other.sum_delta2_noise_;
    sum_noise2_ += other.sum_noise2_;
}


void BinAccumulator::subtract(
    const BinAccumulator& other
) noexcept
{
    count_ -= other.count_;
    sum_delta_squared_ -= other.sum_delta_squared_;
    sum_abs_delta_ -= other.sum_abs_delta_;
    sum_noise_ -= other.sum_noise_;
    sum_delta4_ -= other.sum_delta4_;
    sum_delta2_noise_ -= other.sum_delta2_noise_;
    sum_noise2_ -= other.sum_noise2_;
}


std::size_t
BinAccumulator::count() const noexcept
{
    return count_;
}


double
BinAccumulator::sum_delta_squared() const noexcept
{
    return sum_delta_squared_;
}


double
BinAccumulator::sum_abs_delta() const noexcept
{
    return sum_abs_delta_;
}


double
BinAccumulator::sum_noise() const noexcept
{
    return sum_noise_;
}


double
BinAccumulator::sum_delta4() const noexcept
{
    return sum_delta4_;
}


double
BinAccumulator::sum_delta2_noise() const noexcept
{
    return sum_delta2_noise_;
}


double
BinAccumulator::sum_noise2() const noexcept
{
    return sum_noise2_;
}


double
BinAccumulator::sf_squared(
    SFMethod method
) const noexcept
{
    if (count_ == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double count =
        static_cast<double>(count_);

    switch (method) {
        case SFMethod::SecondOrder:
            // <delta^2> - <noise>
            return (
                sum_delta_squared_ -
                sum_noise_
            ) / count;

        case SFMethod::SecondOrderNoNoise:
            // <delta^2>
            return sum_delta_squared_ / count;

        case SFMethod::MeanAbsoluteDeviation: {
            // pi/2 * <|delta|>^2 - <noise>
            const double mean_abs_delta =
                sum_abs_delta_ / count;

            return (
                kPi / 2.0 *
                mean_abs_delta *
                mean_abs_delta
            ) - (
                sum_noise_ / count
            );
        }

        case SFMethod::MeanAbsoluteDeviationNoNoise: {
            // pi/2 * <|delta|>^2
            const double mean_abs_delta =
                sum_abs_delta_ / count;

            return (
                kPi / 2.0 *
                mean_abs_delta *
                mean_abs_delta
            );
        }
    }

    // Unreachable; keeps the compiler quiet about missing returns.
    return std::numeric_limits<double>::quiet_NaN();
}


double
BinAccumulator::sf(
    SFMethod method
) const noexcept
{
    const double value = sf_squared(method);

    // SF^2 < 0 -> SF = NaN
    if (!std::isfinite(value) || value < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return std::sqrt(value);
}

} // namespace esf
} // namespace agnsf
