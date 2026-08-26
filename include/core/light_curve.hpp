#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace agnsf {

class LightCurveView;

/**
 * Owning time-series data shared by all analysis modules.
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

    /**
     * Non-owning view over this light curve's data.
     *
     * The returned view is valid only while this LightCurve is
     * alive and unmodified.
     */
    LightCurveView view() const noexcept;

private:
    std::vector<double> time_;
    std::vector<double> value_;
    std::vector<double> error_;
};


/**
 * Non-owning view over the arrays of a LightCurve.
 *
 * LightCurveView stores raw pointers only; it does not own or copy
 * any data. The pointed-to buffers must outlive the view. Create a
 * view with LightCurve::view(), or directly from external buffers
 * (e.g. NumPy arrays) whose lifetime is managed by the caller.
 */
class LightCurveView {
public:
    LightCurveView(
        const double* time,
        const double* value,
        const double* error,
        std::size_t size
    ) noexcept;

    std::size_t size() const noexcept;

    const double* time_data() const noexcept;
    const double* value_data() const noexcept;
    const double* error_data() const noexcept;

    std::uintptr_t time_address() const noexcept;
    std::uintptr_t value_address() const noexcept;
    std::uintptr_t error_address() const noexcept;

private:
    const double* time_;
    const double* value_;
    const double* error_;
    std::size_t size_;
};


inline LightCurveView LightCurve::view() const noexcept
{
    return LightCurveView(
        time_data(),
        value_data(),
        error_data(),
        size()
    );
}

} // namespace agnsf
