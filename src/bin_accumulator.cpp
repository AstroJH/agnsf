#include <cmath>
#include <limits>
#include <esf/bin_accumulator.hpp>
#include <iostream>


namespace esf {
    
void BinAccumulator::add(
    double delta,
    double error_i,
    double error_j
) noexcept
{
    const double delta_squared =
        delta * delta;

    const double noise =
        error_i * error_i +
        error_j * error_j;

    sum_delta_squared_ +=
        delta_squared;

    sum_noise_ +=
        noise;

    ++count_;
}


void BinAccumulator::merge(
    const BinAccumulator& other
) noexcept
{
    count_ += other.count_;

    sum_delta_squared_ +=
        other.sum_delta_squared_;

    sum_noise_ +=
        other.sum_noise_;
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
BinAccumulator::sum_noise() const noexcept
{
    return sum_noise_;
}


double
BinAccumulator::sf_squared() const noexcept
{
    if (count_ == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return (
        sum_delta_squared_ -
        sum_noise_
    ) / static_cast<double>(count_);
}


double
BinAccumulator::sf() const noexcept
{
    const double value =
        sf_squared();

    if (!std::isfinite(value) || value < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return std::sqrt(value);
}

} // namespace esf