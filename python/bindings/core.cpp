#include "common.hpp"
#include <esf/lag_bins.hpp>
#include <esf/sf_result.hpp>

void bind_core(py::module_& m)
{
    py::class_<agnsf::esf::SFBinResult>(m, "SFBinResult")
            .def_readonly(
                "count",
                &agnsf::esf::SFBinResult::count
            )
            .def_readonly(
                "sf_squared",
                &agnsf::esf::SFBinResult::sf_squared
            )
            .def_readonly(
                "sf",
                &agnsf::esf::SFBinResult::sf
            )
            .def_readonly(
                "measurement",
                &agnsf::esf::SFBinResult::measurement
            )
            .def_readonly(
                "within",
                &agnsf::esf::SFBinResult::within
            )
            .def_readonly(
                "sampling",
                &agnsf::esf::SFBinResult::sampling
            );
    
    
        // ------------------------------------------------------------------
        // SFResult
        // ------------------------------------------------------------------
    py::class_<agnsf::esf::SFResult>(m, "SFResult")
            .def(
                "__len__",
                &agnsf::esf::SFResult::size
            )
            .def(
                "__getitem__",
                [](const agnsf::esf::SFResult& result, py::ssize_t index)
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
                &agnsf::esf::SFResult::bins
            );
    
    
        // ------------------------------------------------------------------
        // LagBins
        // ------------------------------------------------------------------
    py::class_<agnsf::esf::LagBins> lag_bins(m, "LagBins");
    py::enum_<agnsf::esf::LagBins::GridType>(
            lag_bins,
            "GridType"
        )
            .value(
                "Custom",
                agnsf::esf::LagBins::GridType::Custom
            )
            .value(
                "Linear",
                agnsf::esf::LagBins::GridType::Linear
            )
            .value(
                "Logarithmic",
                agnsf::esf::LagBins::GridType::Logarithmic
            );
    
        lag_bins
            .def(
                py::init<std::vector<double>>(),
                py::arg("edges")
            )
            .def_static(
                "linear",
                &agnsf::esf::LagBins::linear,
                py::arg("min"),
                py::arg("max"),
                py::arg("step")
            )
            .def_static(
                "logarithmic",
                &agnsf::esf::LagBins::logarithmic,
                py::arg("min"),
                py::arg("max"),
                py::arg("step")
            )
            .def(
                "__len__",
                &agnsf::esf::LagBins::size
            )
            .def(
                "__getitem__",
                [](const agnsf::esf::LagBins& bins, py::ssize_t index)
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
                &agnsf::esf::LagBins::grid_type
            )
            .def_property_readonly(
                "min",
                &agnsf::esf::LagBins::min
            )
            .def_property_readonly(
                "max",
                &agnsf::esf::LagBins::max
            )
            .def_property_readonly(
                "size",
                &agnsf::esf::LagBins::size
            )
            .def_property_readonly(
                "edges",
                &agnsf::esf::LagBins::edges
            )
            .def(
                "contains",
                &agnsf::esf::LagBins::contains,
                py::arg("lag")
            )
            .def(
                "index",
                &agnsf::esf::LagBins::index,
                py::arg("lag")
            );
    
    
        // ------------------------------------------------------------------
        // LightCurve
        // ------------------------------------------------------------------
    py::class_<agnsf::LightCurve>(m, "LightCurve")
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
                &agnsf::LightCurve::size
            )
            .def_property_readonly(
                "time",
                &agnsf::LightCurve::time
            )
            .def_property_readonly(
                "value",
                &agnsf::LightCurve::value
            )
            .def_property_readonly(
                "error",
                &agnsf::LightCurve::error
            );
    
    
        // ------------------------------------------------------------------
        // LightCurveView
        // ------------------------------------------------------------------
    py::class_<agnsf::LightCurveView>(m, "LightCurveView")
            .def_property_readonly(
                "size",
                &agnsf::LightCurveView::size
            )
            .def_property_readonly(
                "time_address",
                &agnsf::LightCurveView::time_address
            )
            .def_property_readonly(
                "value_address",
                &agnsf::LightCurveView::value_address
            )
            .def_property_readonly(
                "error_address",
                &agnsf::LightCurveView::error_address
            );
    
    
        // ------------------------------------------------------------------
        // SFMethod
        // ------------------------------------------------------------------
}
