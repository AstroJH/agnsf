#pragma once

#include <cstddef>
#include <cstdint>

namespace esf {

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

} // namespace esf