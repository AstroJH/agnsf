#include <esf/light_curve_view.hpp>

namespace esf {

LightCurveView::LightCurveView(
    const double* time,
    const double* value,
    const double* error,
    std::size_t size
) noexcept
    : time_(time),
      value_(value),
      error_(error),
      size_(size)
{
}

std::size_t LightCurveView::size() const noexcept
{
    return size_;
}

const double* LightCurveView::time_data() const noexcept
{
    return time_;
}

const double* LightCurveView::value_data() const noexcept
{
    return value_;
}

const double* LightCurveView::error_data() const noexcept
{
    return error_;
}

std::uintptr_t LightCurveView::time_address() const noexcept
{
    return reinterpret_cast<std::uintptr_t>(time_);
}

std::uintptr_t LightCurveView::value_address() const noexcept
{
    return reinterpret_cast<std::uintptr_t>(value_);
}

std::uintptr_t LightCurveView::error_address() const noexcept
{
    return reinterpret_cast<std::uintptr_t>(error_);
}

} // namespace esf