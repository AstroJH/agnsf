#include <cstddef>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <esf/lag_bins.hpp>
#include <esf/light_curve.hpp>
#include <esf/sf_calculator.hpp>
#include <esf/sf_result.hpp>

namespace py = pybind11;

PYBIND11_MODULE(_agnsf, m)
{
    m.doc() = "Astronomical structure function analysis";

    // ------------------------------------------------------------------
    // SFBinResult
    // ------------------------------------------------------------------

    py::class_<esf::SFBinResult>(
        m,
        "SFBinResult"
    )
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

    py::class_<esf::SFResult>(
        m,
        "SFResult"
    )
        .def(
            "__len__",
            &esf::SFResult::size
        )
        .def(
            "__getitem__",
            [](const esf::SFResult& result,
               std::size_t index)
            {
                if (index >= result.size()) {
                    throw py::index_error();
                }

                return result.bin(index);
            }
        )
        .def_property_readonly(
            "bins",
            &esf::SFResult::bins
        );

    // ------------------------------------------------------------------
    // LagBins
    // ------------------------------------------------------------------

    py::class_<esf::LagBins> lag_bins(
        m,
        "LagBins"
    );

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

                if (index < 0) {
                    index += size;
                }

                if (index < 0 || index >= size) {
                    throw py::index_error();
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

    py::class_<esf::LightCurve>(
        m,
        "LightCurve"
    )
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
    // Single-light-curve SF
    // ------------------------------------------------------------------

    m.def(
        "sf",
        [](
            const std::vector<double>& time,
            const std::vector<double>& value,
            const std::vector<double>& error,
            const esf::LagBins& bins
        )
        {
            esf::LightCurve light_curve(
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
        py::arg("bins")
    );
}