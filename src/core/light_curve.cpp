#include <core/light_curve.hpp>

#include <stdexcept>
#include <utility>

namespace agnsf {

LightCurve::LightCurve(
    std::vector<double> time,
    std::vector<double> value,
    std::vector<double> error
)
    : time_(std::move(time)),
      value_(std::move(value)),
      error_(std::move(error))
{
    if (time_.size() != value_.size() ||
        time_.size() != error_.size()) {
        throw std::invalid_argument(
            "LightCurve arrays must have the same size."
        );
    }
}


std::size_t LightCurve::size() const noexcept
{
    return time_.size();
}


const std::vector<double>&
LightCurve::time() const noexcept
{
    return time_;
}


const std::vector<double>&
LightCurve::value() const noexcept
{
    return value_;
}


const std::vector<double>&
LightCurve::error() const noexcept
{
    return error_;
}


const double* LightCurve::time_data() const noexcept
{
    return time_.data();
}


const double* LightCurve::value_data() const noexcept
{
    return value_.data();
}


const double* LightCurve::error_data() const noexcept
{
    return error_.data();
}


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

} // namespace agnsf
