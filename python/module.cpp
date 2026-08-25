#include <cstddef>
#include <cstdint>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <esf/lag_bins.hpp>
#include <esf/light_curve.hpp>
#include <esf/light_curve_view.hpp>
#include <esf/sf_calculator.hpp>
#include <esf/pooled_esf_calculator.hpp>
#include <esf/sf_result.hpp>
#include <esf/sf_ensemble_calculator.hpp>

namespace py = pybind11;

namespace {

using Float64Array =
    py::array_t<double, py::array::c_style>;


// ----------------------------------------------------------------------
// Validate three input arrays and construct a zero-copy LightCurveView.
// ----------------------------------------------------------------------

esf::LightCurveView make_light_curve_view(
    const Float64Array& time,
    const Float64Array& value,
    const Float64Array& error
)
{
    const auto time_buffer = time.request();
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

    return esf::LightCurveView(
        time_ptr,
        value_ptr,
        error_ptr,
        static_cast<std::size_t>(time_buffer.size)
    );
}

// ----------------------------------------------------------------------
// Convert a list of light curves into zero-copy LightCurveViews.
//
// LightCurveView only stores raw pointers.  If ensure() must build a
// new float64 C-contiguous array (e.g. from a Python list, a float32
// array, or a non-contiguous slice), that temporary ndarray would be
// destroyed at the end of the loop iteration.  The owned_* vectors
// keep every array alive until the caller finishes using the views.
// ----------------------------------------------------------------------

struct LightCurveBatch {
    std::vector<esf::LightCurveView> views;

    std::vector<py::array_t<double, py::array::c_style>> owned_times;
    std::vector<py::array_t<double, py::array::c_style>> owned_values;
    std::vector<py::array_t<double, py::array::c_style>> owned_errors;
};


LightCurveBatch make_light_curve_batch(
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

        auto time  = py::array_t<double, py::array::c_style>::ensure(times[i]);
        auto value = py::array_t<double, py::array::c_style>::ensure(values[i]);
        auto error = py::array_t<double, py::array::c_style>::ensure(errors[i]);

        batch.owned_times.emplace_back(time);
        batch.owned_values.emplace_back(value);
        batch.owned_errors.emplace_back(error);

        if (!time || !value || !error) {
            throw py::type_error(
                "all light curves must be convertible "
                "to C-contiguous float64 arrays"
            );
        }

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

        batch.views.emplace_back(
            static_cast<const double*>(time_buffer.ptr),
            static_cast<const double*>(value_buffer.ptr),
            static_cast<const double*>(error_buffer.ptr),
            static_cast<std::size_t>(time_buffer.size)
        );
    }

    return batch;
}

} // namespace


PYBIND11_MODULE(_agnsf, m)
{
    m.doc() = "Astronomical structure function analysis";


    // ------------------------------------------------------------------
    // SFBinResult
    // ------------------------------------------------------------------

    py::class_<esf::SFBinResult>(m, "SFBinResult")
        .def_readonly(
            "count",
            &esf::SFBinResult::count
        )
        .def_readonly(
            "sf_squared",
            &esf::SFBinResult::sf_squared
        )
        .def_readonly(
            "sf",
            &esf::SFBinResult::sf
        );


    // ------------------------------------------------------------------
    // SFResult
    // ------------------------------------------------------------------

    py::class_<esf::SFResult>(m, "SFResult")
        .def(
            "__len__",
            &esf::SFResult::size
        )
        .def(
            "__getitem__",
            [](const esf::SFResult& result, py::ssize_t index)
            {
                const py::ssize_t size =
                    static_cast<py::ssize_t>(result.size());

                // Python-style negative indexing.
                if (index < 0) {
                    index += size;
                }

                if (index < 0 || index >= size) {
                    throw py::index_error(
                        "SFResult index out of range"
                    );
                }

                return result.bin(
                    static_cast<std::size_t>(index)
                );
            }
        )
        .def_property_readonly(
            "bins",
            &esf::SFResult::bins
        );


    // ------------------------------------------------------------------
    // LagBins
    // ------------------------------------------------------------------

    py::class_<esf::LagBins> lag_bins(m, "LagBins");

    py::enum_<esf::LagBins::GridType>(
        lag_bins,
        "GridType"
    )
        .value(
            "Custom",
            esf::LagBins::GridType::Custom
        )
        .value(
            "Linear",
            esf::LagBins::GridType::Linear
        )
        .value(
            "Logarithmic",
            esf::LagBins::GridType::Logarithmic
        );

    lag_bins
        .def(
            py::init<std::vector<double>>(),
            py::arg("edges")
        )
        .def_static(
            "linear",
            &esf::LagBins::linear,
            py::arg("min"),
            py::arg("max"),
            py::arg("step")
        )
        .def_static(
            "logarithmic",
            &esf::LagBins::logarithmic,
            py::arg("min"),
            py::arg("max"),
            py::arg("step")
        )
        .def(
            "__len__",
            &esf::LagBins::size
        )
        .def(
            "__getitem__",
            [](const esf::LagBins& bins, py::ssize_t index)
            {
                const py::ssize_t size =
                    static_cast<py::ssize_t>(bins.size());

                // Python-style negative indexing.
                if (index < 0) {
                    index += size;
                }

                if (index < 0 || index >= size) {
                    throw py::index_error(
                        "LagBins index out of range"
                    );
                }

                const auto& edges = bins.edges();

                return py::make_tuple(
                    edges[static_cast<std::size_t>(index)],
                    edges[static_cast<std::size_t>(index) + 1]
                );
            }
        )
        .def_property_readonly(
            "grid_type",
            &esf::LagBins::grid_type
        )
        .def_property_readonly(
            "min",
            &esf::LagBins::min
        )
        .def_property_readonly(
            "max",
            &esf::LagBins::max
        )
        .def_property_readonly(
            "size",
            &esf::LagBins::size
        )
        .def_property_readonly(
            "edges",
            &esf::LagBins::edges
        )
        .def(
            "contains",
            &esf::LagBins::contains,
            py::arg("lag")
        )
        .def(
            "index",
            &esf::LagBins::index,
            py::arg("lag")
        );


    // ------------------------------------------------------------------
    // LightCurve
    // ------------------------------------------------------------------

    py::class_<esf::LightCurve>(m, "LightCurve")
        .def(
            py::init<
                std::vector<double>,
                std::vector<double>,
                std::vector<double>
            >(),
            py::arg("time"),
            py::arg("value"),
            py::arg("error")
        )
        .def_property_readonly(
            "size",
            &esf::LightCurve::size
        )
        .def_property_readonly(
            "time",
            &esf::LightCurve::time
        )
        .def_property_readonly(
            "value",
            &esf::LightCurve::value
        )
        .def_property_readonly(
            "error",
            &esf::LightCurve::error
        );


    // ------------------------------------------------------------------
    // LightCurveView
    // ------------------------------------------------------------------

    py::class_<esf::LightCurveView>(m, "LightCurveView")
        .def_property_readonly(
            "size",
            &esf::LightCurveView::size
        )
        .def_property_readonly(
            "time_address",
            &esf::LightCurveView::time_address
        )
        .def_property_readonly(
            "value_address",
            &esf::LightCurveView::value_address
        )
        .def_property_readonly(
            "error_address",
            &esf::LightCurveView::error_address
        );


    // ------------------------------------------------------------------
    // Single-light-curve SF
    //
    // NumPy arrays are accepted only when they are:
    //
    //   - float64
    //   - one-dimensional
    //   - C-contiguous
    //
    // Therefore no implicit dtype conversion or contiguous-copy is
    // performed here.
    // ------------------------------------------------------------------

    m.def(
        "sf",
        [](
            const Float64Array& time,
            const Float64Array& value,
            const Float64Array& error,
            const esf::LagBins& bins
        )
        {
            const auto light_curve =
                make_light_curve_view(
                    time,
                    value,
                    error
                );

            esf::SFCalculator calculator;

            return calculator.calculate(
                light_curve,
                bins
            );
        },
        py::arg("time"),
        py::arg("value"),
        py::arg("error"),
        py::arg("bins"),
        R"pbdoc(
Calculate the structure function of a single light curve.

The input arrays must be one-dimensional, float64, and
C-contiguous. No data copy is performed by this interface.
        )pbdoc"
    );

    m.def(
        "pooled_sf",
        [](
            const py::list& times,
            const py::list& values,
            const py::list& errors,
            const esf::LagBins& bins
        )
        {
            const auto batch =
                make_light_curve_batch(
                    times,
                    values,
                    errors
                );

            esf::PooledESFCalculator calculator;

            return calculator.calculate(
                batch.views,
                bins
            );
        },
        py::arg("times"),
        py::arg("values"),
        py::arg("errors"),
        py::arg("bins")
    );


    // ------------------------------------------------------------------
    // EnsembleMethod
    // ------------------------------------------------------------------

    py::enum_<esf::SFEnsembleCalculator::Method>(
        m,
        "EnsembleMethod"
    )
        .value(
            "SqrtMeanSquared",
            esf::SFEnsembleCalculator::Method::SqrtMeanSquared
        )
        .value(
            "MeanSf",
            esf::SFEnsembleCalculator::Method::MeanSf
        );


    // ------------------------------------------------------------------
    // Aggregated ensemble SF
    //
    // Computes an SF for each light curve and combines the individual
    // SF values according to the chosen method:
    //
    //   SqrtMeanSquared:  ESF = sqrt( <SF^2> )
    //   MeanSf:           ESF = <SF>
    //
    // The same zero-copy / lifetime-safe view handling as pooled_sf
    // is used.
    // ------------------------------------------------------------------

    m.def(
        "ensemble_sf",
        [](
            const py::list& times,
            const py::list& values,
            const py::list& errors,
            const esf::LagBins& bins,
            esf::SFEnsembleCalculator::Method method
        )
        {
            const auto batch =
                make_light_curve_batch(
                    times,
                    values,
                    errors
                );

            esf::SFEnsembleCalculator calculator;

            return calculator.calculate(
                batch.views,
                bins,
                method
            );
        },
        py::arg("times"),
        py::arg("values"),
        py::arg("errors"),
        py::arg("bins"),
        py::arg("method") =
            esf::SFEnsembleCalculator::Method::SqrtMeanSquared,
        R"pbdoc(
Calculate the aggregated ensemble structure function.

An SF is computed independently for each light curve and the
individual SF values are combined per lag bin.

method: EnsembleMethod.SqrtMeanSquared (default)
            ESF(tau) = sqrt( <SF_k^2(tau)>_k )
        EnsembleMethod.MeanSf
            ESF(tau) = <SF_k(tau)>_k

Only light curves with a finite contribution are included in
each lag bin. Each contributing light curve is weighted equally.
        )pbdoc"
    );


    // ------------------------------------------------------------------
    // Inspect NumPy inputs
    //
    // This function is intentionally separate from sf().
    // It accepts arbitrary Python objects and reports whether pybind11
    // would need to construct a new float64 C-contiguous array.
    // ------------------------------------------------------------------

    m.def(
        "inspect",
        [](
            py::object time,
            py::object value,
            py::object error
        )
        {
            auto time_array =
                py::array_t<
                    double,
                    py::array::c_style |
                    py::array::forcecast
                >::ensure(time);

            auto value_array =
                py::array_t<
                    double,
                    py::array::c_style |
                    py::array::forcecast
                >::ensure(value);

            auto error_array =
                py::array_t<
                    double,
                    py::array::c_style |
                    py::array::forcecast
                >::ensure(error);

            if (!time_array ||
                !value_array ||
                !error_array) {

                throw py::type_error(
                    "inputs cannot be converted to float64 "
                    "C-contiguous arrays"
                );
            }

            py::dict result;

            result["time_address"] =
                reinterpret_cast<std::uintptr_t>(
                    time_array.data()
                );

            result["value_address"] =
                reinterpret_cast<std::uintptr_t>(
                    value_array.data()
                );

            result["error_address"] =
                reinterpret_cast<std::uintptr_t>(
                    error_array.data()
                );

            result["time_copied"] =
                !time.is(time_array);

            result["value_copied"] =
                !value.is(value_array);

            result["error_copied"] =
                !error.is(error_array);

            return result;
        },
        py::arg("time"),
        py::arg("value"),
        py::arg("error"),
        R"pbdoc(
Inspect whether the supplied inputs require conversion to
float64 C-contiguous NumPy arrays.

This function may create temporary arrays when conversion is
necessary. It is intended for diagnostics and does not perform
the structure-function calculation.
        )pbdoc"
    );
}