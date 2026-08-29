#pragma once

#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <core/light_curve.hpp>
#include <core/uncertainty.hpp>

namespace py = pybind11;

namespace agnsf::python {

using Float64Array =
    py::array_t<double, py::array::c_style>;

// Ensure that a Python object can be represented as a C-contiguous
// float64 NumPy array.
inline Float64Array ensure_float64_array(
    const py::handle& object
)
{
    auto array = Float64Array::ensure(object);

    if (!array) {
        throw py::type_error(
            "input must be convertible to "
            "a C-contiguous float64 array"
        );
    }

    return array;
}


// Construct an owning LightCurve by copying the input arrays.
inline agnsf::LightCurve make_light_curve(
    const Float64Array& time,
    const Float64Array& value,
    const Float64Array& error
)
{
    return agnsf::LightCurve(
        py::cast<std::vector<double>>(time),
        py::cast<std::vector<double>>(value),
        py::cast<std::vector<double>>(error)
    );
}

// Construct a zero-copy LightCurveView.
//
// LightCurveView does not own the underlying data; the input arrays
// must remain alive for the lifetime of the view.
inline agnsf::LightCurveView make_light_curve_view(
    const Float64Array& time,
    const Float64Array& value,
    const Float64Array& error
)
{
    const auto time_buffer  = time.request();
    const auto value_buffer = value.request();
    const auto error_buffer = error.request();

    if (time_buffer.ndim  != 1 ||
        value_buffer.ndim != 1 ||
        error_buffer.ndim != 1) {

        throw py::value_error(
            "time, value, and error must be 1-dimensional"
        );
    }

    if (time_buffer.size != value_buffer.size ||
        time_buffer.size != error_buffer.size) {

        throw py::value_error(
            "time, value, and error must have the same length"
        );
    }

    const auto* time_ptr =
        static_cast<const double*>(time_buffer.ptr);

    const auto* value_ptr =
        static_cast<const double*>(value_buffer.ptr);

    const auto* error_ptr =
        static_cast<const double*>(error_buffer.ptr);

    return agnsf::LightCurveView(
        time_ptr,
        value_ptr,
        error_ptr,
        static_cast<std::size_t>(time_buffer.size)
    );
}

struct LightCurveBatch {
    std::vector<agnsf::LightCurveView> views;
    
    // Keep converted NumPy arrays alive for the lifetime of the views.
    std::vector<Float64Array> owned_times;
    std::vector<Float64Array> owned_values;
    std::vector<Float64Array> owned_errors;
};


// Convert a list of light curves into zero-copy LightCurveViews.
//
// LightCurveView only stores raw pointers.  If ensure() must build a
// new float64 C-contiguous array (e.g. from a Python list, a float32
// array, or a non-contiguous slice), that temporary ndarray would be
// destroyed at the end of the loop iteration.  The owned_* vectors
// keep every array alive until the caller finishes using the views.
inline LightCurveBatch make_light_curve_batch(
    const py::list& times,
    const py::list& values,
    const py::list& errors
)
{
    if (times.size() != values.size() ||
        times.size() != errors.size()) {

        throw py::value_error(
            "times, values, and errors must contain "
            "the same number of light curves"
        );
    }

    LightCurveBatch batch;
    batch.views.reserve(times.size());
    batch.owned_times.reserve(times.size());
    batch.owned_values.reserve(values.size());
    batch.owned_errors.reserve(errors.size());

    for (std::size_t i = 0; i < times.size(); ++i) {
        // >> OWNERSHIP WARNING <<
        // If ensure() created a copy, the temporary ndarray is destroyed
        // when this loop iteration ends; LightCurveView only keeps raw pointers.
        // Keep the resulting arrays in owned_* so that their buffers remain valid
        // until calculate() finishes.

        auto time  = ensure_float64_array(times[i]);
        auto value = ensure_float64_array(values[i]);
        auto error = ensure_float64_array(errors[i]);

        auto time_buffer  = time.request();
        auto value_buffer = value.request();
        auto error_buffer = error.request();

        if (time_buffer.ndim  != 1 ||
            value_buffer.ndim != 1 ||
            error_buffer.ndim != 1) {

            throw py::value_error(
                "time, value, and error must be "
                "1-dimensional"
            );
        }

        if (time_buffer.size != value_buffer.size ||
            time_buffer.size != error_buffer.size) {

            throw py::value_error(
                "time, value, and error must have "
                "the same length for each light curve"
            );
        }

        batch.owned_times.emplace_back(time);
        batch.owned_values.emplace_back(value);
        batch.owned_errors.emplace_back(error);

        batch.views.emplace_back(
            static_cast<const double*>(time_buffer.ptr),
            static_cast<const double*>(value_buffer.ptr),
            static_cast<const double*>(error_buffer.ptr),
            static_cast<std::size_t>(time_buffer.size)
        );
    }

    return batch;
}

} // namespace agnsf::python


// Registration entry points for each binding group.
void bind_core(py::module_& m);
void bind_uncertainty(py::module_& m);
void bind_sf(py::module_& m);
void bind_esf(py::module_& m);
void bind_io(py::module_& m);
void bind_timedelay(py::module_& m);
