#pragma once

#include <cstddef>
#include <vector>

namespace esf {

/**
 * Time-series data used by the SF/ESF calculators.
 *
 * Contract:
 *
 * 1. time, value, and error must have the same length.
 * 2. time must be sorted in non-decreasing order.
 * 3. time, value, and error must contain finite values.
 * 4. This class does not sort the input or validate
 *    the ordering or finiteness of the input data.
 */
class LightCurve {
public:
    LightCurve(
        std::vector<double> time,
        std::vector<double> value,
        std::vector<double> error
    );

    std::size_t size() const noexcept;

    const std::vector<double>& time() const noexcept;
    const std::vector<double>& value() const noexcept;
    const std::vector<double>& error() const noexcept;

    const double* time_data() const noexcept;
    const double* value_data() const noexcept;
    const double* error_data() const noexcept;

private:
    std::vector<double> time_;
    std::vector<double> value_;
    std::vector<double> error_;
};

} // namespace esf