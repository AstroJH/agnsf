#include <esf/light_curve.hpp>

#include <stdexcept>
#include <utility>

namespace esf {

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

} // namespace esf