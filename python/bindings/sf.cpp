#include "common.hpp"
#include <esf/sf_method.hpp>
#include <esf/sf_calculator.hpp>

using namespace agnsf::python;

void bind_sf(py::module_& m)
{
    // ------------------------------------------------------------------
    // SFMethod
    // ------------------------------------------------------------------
    py::enum_<agnsf::esf::SFMethod>(
            m,
            "SFMethod"
        )
            .value(
                "SecondOrder",
                agnsf::esf::SFMethod::SecondOrder
            )
            .value(
                "SecondOrderNoNoise",
                agnsf::esf::SFMethod::SecondOrderNoNoise
            )
            .value(
                "MeanAbsoluteDeviation",
                agnsf::esf::SFMethod::MeanAbsoluteDeviation
            )
            .value(
                "MeanAbsoluteDeviationNoNoise",
                agnsf::esf::SFMethod::MeanAbsoluteDeviationNoNoise
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
                const agnsf::esf::LagBins& bins,
                agnsf::esf::SFMethod method,
                const agnsf::UncertaintyConfig& uncertainty,
                double redshift
            )
            {
                const auto light_curve =
                    make_light_curve_view(
                        time,
                        value,
                        error
                    );
    
                agnsf::esf::SFCalculator calculator;
    
                return calculator.calculate(
                    light_curve,
                    bins,
                    method,
                    uncertainty,
                    redshift
                );
            },
            py::arg("time"),
            py::arg("value"),
            py::arg("error"),
            py::arg("bins"),
            py::arg("method") =
                agnsf::esf::SFMethod::SecondOrder,
            py::arg("uncertainty") =
                agnsf::UncertaintyConfig(),
            py::arg("redshift") = 0.0,
            R"pbdoc(
    Calculate the structure function of a single light curve.
    
    The input arrays must be one-dimensional, float64, and
    C-contiguous. No data copy is performed by this interface.
    
    method: SFMethod.SecondOrder (default)
                SF^2 = <delta^2> - <sigma_i^2 + sigma_j^2>
            SFMethod.SecondOrderNoNoise
                SF^2 = <delta^2>
            SFMethod.MeanAbsoluteDeviation
                SF^2 = pi/2 * <|delta|>^2 - <sigma_i^2 + sigma_j^2>
            SFMethod.MeanAbsoluteDeviationNoNoise
                SF^2 = pi/2 * <|delta|>^2
    
    uncertainty: UncertaintyConfig; measurement=Analytic estimates the
                within-bin standard error of the mean propagated to sf.
                Sampling is not defined for a single light curve.
    
    redshift: source redshift z; lags are converted to the rest frame
                (dt_rest = dt_obs / (1 + z)).
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
