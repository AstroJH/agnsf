#include <cstddef>
#include <cstdint>
#include <sstream>
#include <string>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <esf/lag_bins.hpp>
#include <core/light_curve.hpp>
#include <core/light_curve.hpp>
#include <esf/sf_calculator.hpp>
#include <esf/pooled_esf_calculator.hpp>
#include <esf/sf_result.hpp>
#include <esf/sf_ensemble_calculator.hpp>
#include <esf/sf_io.hpp>
#include <esf/sf_uncertainty.hpp>

#include <io/light_curve.hpp>
#include <io/path_list.hpp>
#include <io/table_writer.hpp>

namespace py = pybind11;

namespace {

using Float64Array =
    py::array_t<double, py::array::c_style>;


// ----------------------------------------------------------------------
// Validate three input arrays and construct a zero-copy LightCurveView.
// ----------------------------------------------------------------------

agnsf::LightCurveView make_light_curve_view(
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

    return agnsf::LightCurveView(
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
    std::vector<agnsf::LightCurveView> views;

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
    // SFUncertainty
    // ------------------------------------------------------------------

    py::class_<agnsf::esf::SFUncertainty>(m, "SFUncertainty")
        .def_readonly(
            "lower",
            &agnsf::esf::SFUncertainty::lower
        )
        .def_readonly(
            "upper",
            &agnsf::esf::SFUncertainty::upper
        )
        .def_property_readonly(
            "estimated",
            &agnsf::esf::SFUncertainty::estimated
        )
        .def(
            "__repr__",
            [](const agnsf::esf::SFUncertainty& uncertainty)
            {
                std::ostringstream out;
                out << "SFUncertainty(lower="
                    << uncertainty.lower
                    << ", upper="
                    << uncertainty.upper
                    << ")";
                return out.str();
            }
        );


    // ------------------------------------------------------------------
    // UncertaintyMethod
    // ------------------------------------------------------------------

    py::enum_<agnsf::esf::UncertaintyMethod>(
        m,
        "UncertaintyMethod"
    )
        .value(
            "Off",
            agnsf::esf::UncertaintyMethod::Off
        )
        .value(
            "Analytic",
            agnsf::esf::UncertaintyMethod::Analytic
        )
        .value(
            "Jackknife",
            agnsf::esf::UncertaintyMethod::Jackknife
        )
        .value(
            "Bootstrap",
            agnsf::esf::UncertaintyMethod::Bootstrap
        );


    // ------------------------------------------------------------------
    // UncertaintyConfig
    // ------------------------------------------------------------------

    py::class_<agnsf::esf::UncertaintyConfig>(
        m,
        "UncertaintyConfig"
    )
        .def(
            py::init(
                [](
                    agnsf::esf::UncertaintyMethod measurement,
                    agnsf::esf::UncertaintyMethod sampling,
                    std::size_t n_bootstrap,
                    std::uint32_t bootstrap_seed
                )
                {
                    agnsf::esf::UncertaintyConfig config;

                    config.measurement = measurement;
                    config.sampling = sampling;
                    config.n_bootstrap = n_bootstrap;
                    config.bootstrap_seed = bootstrap_seed;

                    return config;
                }
            ),
            py::arg("measurement") =
                agnsf::esf::UncertaintyMethod::Off,
            py::arg("sampling") =
                agnsf::esf::UncertaintyMethod::Off,
            py::arg("n_bootstrap") =
                std::size_t{100},
            py::arg("bootstrap_seed") =
                std::uint32_t{0}
        )
        .def_readwrite(
            "measurement",
            &agnsf::esf::UncertaintyConfig::measurement
        )
        .def_readwrite(
            "sampling",
            &agnsf::esf::UncertaintyConfig::sampling
        )
        .def_readwrite(
            "n_bootstrap",
            &agnsf::esf::UncertaintyConfig::n_bootstrap
        )
        .def_readwrite(
            "bootstrap_seed",
            &agnsf::esf::UncertaintyConfig::bootstrap_seed
        )
        .def(
            "__repr__",
            [](const agnsf::esf::UncertaintyConfig& config)
            {
                const auto name =
                    [](agnsf::esf::UncertaintyMethod method)
                    {
                        switch (method) {
                            case agnsf::esf::UncertaintyMethod::Off:
                                return "Off";
                            case agnsf::esf::UncertaintyMethod::Analytic:
                                return "Analytic";
                            case agnsf::esf::UncertaintyMethod::Jackknife:
                                return "Jackknife";
                            case agnsf::esf::UncertaintyMethod::Bootstrap:
                                return "Bootstrap";
                        }
                        return "?";
                    };

                std::ostringstream out;
                out << "UncertaintyConfig(measurement="
                    << name(config.measurement)
                    << ", sampling="
                    << name(config.sampling)
                    << ", n_bootstrap="
                    << config.n_bootstrap
                    << ", bootstrap_seed="
                    << config.bootstrap_seed
                    << ")";
                return out.str();
            }
        );


    // ------------------------------------------------------------------
    // SFBinResult
    // ------------------------------------------------------------------

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
            const agnsf::esf::UncertaintyConfig& uncertainty,
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
            agnsf::esf::UncertaintyConfig(),
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

    m.def(
        "pooled_sf",
        [](
            const py::list& times,
            const py::list& values,
            const py::list& errors,
            const agnsf::esf::LagBins& bins,
            agnsf::esf::SFMethod method,
            const agnsf::esf::UncertaintyConfig& uncertainty,
            const py::object& redshift
        )
        {
            const auto batch =
                make_light_curve_batch(
                    times,
                    values,
                    errors
                );

            agnsf::esf::PooledESFCalculator calculator;

            if (py::isinstance<py::float_>(redshift) ||
                py::isinstance<py::int_>(redshift)) {

                // Scalar: the same z for every light curve.
                return calculator.calculate(
                    batch.views,
                    bins,
                    method,
                    uncertainty,
                    redshift.cast<double>()
                );
            }

            // Sequence: one redshift per light curve.
            const std::vector<double> redshifts =
                redshift.cast<std::vector<double>>();

            return calculator.calculate(
                batch.views,
                bins,
                method,
                uncertainty,
                redshifts
            );
        },
        py::arg("times"),
        py::arg("values"),
        py::arg("errors"),
        py::arg("bins"),
        py::arg("method") =
            agnsf::esf::SFMethod::SecondOrder,
        py::arg("uncertainty") =
            agnsf::esf::UncertaintyConfig(),
        py::arg("redshift") = 0.0,
        R"pbdoc(
Calculate the pooled ensemble structure function.

All valid pairs from all light curves are pooled into the same
lag bins. See sf() for the available SFMethod estimators.

uncertainty: UncertaintyConfig; measurement=Analytic uses the pooled
            pair statistics; sampling supports Jackknife / Bootstrap
            (curve-level). Analytic sampling is not defined here.

redshift: source redshift z — a scalar applies the same z to all
            light curves; a sequence provides one z per light curve.
            Lags are converted to the rest frame
            (dt_rest = dt_obs / (1 + z)).
        )pbdoc"
    );


    // ------------------------------------------------------------------
    // EnsembleMethod
    // ------------------------------------------------------------------

    py::enum_<agnsf::esf::SFEnsembleCalculator::Method>(
        m,
        "EnsembleMethod"
    )
        .value(
            "SqrtMeanSquared",
            agnsf::esf::SFEnsembleCalculator::Method::SqrtMeanSquared
        )
        .value(
            "MeanSf",
            agnsf::esf::SFEnsembleCalculator::Method::MeanSf
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
            const agnsf::esf::LagBins& bins,
            agnsf::esf::SFEnsembleCalculator::Method method,
            agnsf::esf::SFMethod sf_method,
            const agnsf::esf::UncertaintyConfig& uncertainty,
            const py::object& redshift
        )
        {
            const auto batch =
                make_light_curve_batch(
                    times,
                    values,
                    errors
                );

            agnsf::esf::SFEnsembleCalculator calculator;

            if (py::isinstance<py::float_>(redshift) ||
                py::isinstance<py::int_>(redshift)) {

                // Scalar: the same z for every light curve.
                return calculator.calculate(
                    batch.views,
                    bins,
                    sf_method,
                    method,
                    uncertainty,
                    redshift.cast<double>()
                );
            }

            // Sequence: one redshift per light curve.
            const std::vector<double> redshifts =
                redshift.cast<std::vector<double>>();

            return calculator.calculate(
                batch.views,
                bins,
                sf_method,
                method,
                uncertainty,
                redshifts
            );
        },
        py::arg("times"),
        py::arg("values"),
        py::arg("errors"),
        py::arg("bins"),
        py::arg("method") =
            agnsf::esf::SFEnsembleCalculator::Method::SqrtMeanSquared,
        py::arg("sf_method") =
            agnsf::esf::SFMethod::SecondOrder,
        py::arg("uncertainty") =
            agnsf::esf::UncertaintyConfig(),
        py::arg("redshift") = 0.0,
        R"pbdoc(
Calculate the aggregated ensemble structure function.

An SF is computed independently for each light curve using the
chosen SFMethod estimator, and the individual SF values are then
combined per lag bin.

method: EnsembleMethod.SqrtMeanSquared (default)
            ESF(tau) = sqrt( <SF_k^2(tau)>_k )
        EnsembleMethod.MeanSf
            ESF(tau) = <SF_k(tau)>_k

sf_method: SFMethod.SecondOrder (default)
            see sf() for the available estimators.

uncertainty: UncertaintyConfig; measurement=Analytic propagates the
            per-curve measurement; sampling supports Analytic
            (default), Jackknife or Bootstrap.

redshift: source redshift z — a scalar applies the same z to all
            light curves; a sequence provides one z per light curve.
            Lags are converted to the rest frame
            (dt_rest = dt_obs / (1 + z)).

Only light curves with a finite contribution are included in
each lag bin. Each contributing light curve is weighted equally.
        )pbdoc"
    );


    // ------------------------------------------------------------------
    // File IO
    // ------------------------------------------------------------------

    m.def(
        "read_light_curve",
        [](
            const std::string& path,
            const std::string& time,
            const std::string& value,
            const std::string& error
        )
        {
            agnsf::io::ColumnNames columns;
            columns.time = time;
            columns.value = value;
            columns.error = error;

            return agnsf::io::read_light_curve(
                path,
                columns
            );
        },
        py::arg("path"),
        py::arg("time") = "time",
        py::arg("value") = "value",
        py::arg("error") = "error",
        R"pbdoc(
Read a light curve from a CSV or FITS file.

The format is chosen from the file extension (.csv, .fits, .fit,
.fits.gz, .fit.gz). Columns are located by name for files with a
header; otherwise the first three columns are used.

Returns a LightCurve with time, value, and error.
        )pbdoc"
    );

    m.def(
        "read_path_list",
        [](
            const std::string& path
        )
        {
            return agnsf::io::read_path_list(path);
        },
        py::arg("path"),
        R"pbdoc(
Read a list of file paths from a text file (one per line, '#' comments
and blank lines are ignored).
        )pbdoc"
    );

    m.def(
        "read_path_list_with_redshift",
        [](
            const std::string& path
        )
        {
            const auto entries =
                agnsf::io::read_path_list_with_redshift(path);

            py::list result;

            for (const auto& entry : entries) {
                result.append(
                    py::make_tuple(entry.path, entry.redshift)
                );
            }

            return result;
        },
        py::arg("path"),
        R"pbdoc(
Read a list of (path, redshift) entries from a text file.

Each line may have one or two whitespace-separated columns:

  /data/lc1.csv
  /data/lc2.csv 0.5

A missing redshift column means no correction (z = 0). Returns a list
of (path, redshift) tuples.
        )pbdoc"
    );

    m.def(
        "write_table",
        [](
            const std::string& path,
            const std::vector<std::string>& headers,
            const std::vector<std::vector<double>>& columns
        )
        {
            agnsf::io::write_table(path, headers, columns);
        },
        py::arg("path"),
        py::arg("headers"),
        py::arg("columns"),
        R"pbdoc(
Write columns of numbers to a simple text file.

The first line is '# header_1 header_2 ...' followed by one row per
line with space-separated values.
        )pbdoc"
    );


    // ------------------------------------------------------------------
    // File-based SF / ESF convenience functions
    // ------------------------------------------------------------------

    m.def(
        "sf_from_file",
        [](
            const std::string& path,
            const agnsf::esf::LagBins& bins,
            agnsf::esf::SFMethod method,
            const std::string& time,
            const std::string& value,
            const std::string& error,
            const agnsf::esf::UncertaintyConfig& uncertainty,
            double redshift
        )
        {
            agnsf::io::ColumnNames columns;
            columns.time = time;
            columns.value = value;
            columns.error = error;

            return agnsf::esf::sf_from_file(
                path,
                bins,
                method,
                columns,
                uncertainty,
                redshift
            );
        },
        py::arg("path"),
        py::arg("bins"),
        py::arg("method") =
            agnsf::esf::SFMethod::SecondOrder,
        py::arg("time") = "time",
        py::arg("value") = "value",
        py::arg("error") = "error",
        py::arg("uncertainty") =
            agnsf::esf::UncertaintyConfig(),
        py::arg("redshift") = 0.0,
        R"pbdoc(
Calculate the structure function of one light-curve file (CSV or FITS).

path:  light-curve file
bins:  LagBins
method: SFMethod estimator (see sf())
uncertainty: UncertaintyConfig (see sf())
redshift: source redshift z (rest-frame lags)
        )pbdoc"
    );

    m.def(
        "pooled_sf_from_files",
        [](
            const py::list& paths,
            const agnsf::esf::LagBins& bins,
            agnsf::esf::SFMethod method,
            const std::string& time,
            const std::string& value,
            const std::string& error,
            const agnsf::esf::UncertaintyConfig& uncertainty,
            double redshift
        )
        {
            std::vector<std::string> path_vector;
            path_vector.reserve(paths.size());

            for (const auto item : paths) {
                path_vector.push_back(
                    py::cast<std::string>(item)
                );
            }

            agnsf::io::ColumnNames columns;
            columns.time = time;
            columns.value = value;
            columns.error = error;

            return agnsf::esf::pooled_sf_from_files(
                path_vector,
                bins,
                method,
                columns,
                uncertainty,
                redshift
            );
        },
        py::arg("paths"),
        py::arg("bins"),
        py::arg("method") =
            agnsf::esf::SFMethod::SecondOrder,
        py::arg("time") = "time",
        py::arg("value") = "value",
        py::arg("error") = "error",
        py::arg("uncertainty") =
            agnsf::esf::UncertaintyConfig(),
        py::arg("redshift") = 0.0,
        R"pbdoc(
Calculate the pooled ensemble structure function from a list of
light-curve files (CSV or FITS).
        )pbdoc"
    );

    m.def(
        "pooled_sf_from_path_list",
        [](
            const std::string& path_list_file,
            const agnsf::esf::LagBins& bins,
            agnsf::esf::SFMethod method,
            const std::string& time,
            const std::string& value,
            const std::string& error,
            const agnsf::esf::UncertaintyConfig& uncertainty
        )
        {
            agnsf::io::ColumnNames columns;
            columns.time = time;
            columns.value = value;
            columns.error = error;

            return agnsf::esf::pooled_sf_from_path_list(
                path_list_file,
                bins,
                method,
                columns,
                uncertainty
            );
        },
        py::arg("path_list_file"),
        py::arg("bins"),
        py::arg("method") =
            agnsf::esf::SFMethod::SecondOrder,
        py::arg("time") = "time",
        py::arg("value") = "value",
        py::arg("error") = "error",
        py::arg("uncertainty") =
            agnsf::esf::UncertaintyConfig(),
        R"pbdoc(
Calculate the pooled ensemble structure function from a text file that
lists the light-curve paths (one per line).
        )pbdoc"
    );

    m.def(
        "ensemble_sf_from_files",
        [](
            const py::list& paths,
            const agnsf::esf::LagBins& bins,
            agnsf::esf::SFEnsembleCalculator::Method method,
            agnsf::esf::SFMethod sf_method,
            const std::string& time,
            const std::string& value,
            const std::string& error,
            const agnsf::esf::UncertaintyConfig& uncertainty,
            double redshift
        )
        {
            std::vector<std::string> path_vector;
            path_vector.reserve(paths.size());

            for (const auto item : paths) {
                path_vector.push_back(
                    py::cast<std::string>(item)
                );
            }

            agnsf::io::ColumnNames columns;
            columns.time = time;
            columns.value = value;
            columns.error = error;

            return agnsf::esf::ensemble_sf_from_files(
                path_vector,
                bins,
                sf_method,
                method,
                columns,
                uncertainty,
                redshift
            );
        },
        py::arg("paths"),
        py::arg("bins"),
        py::arg("method") =
            agnsf::esf::SFEnsembleCalculator::Method::SqrtMeanSquared,
        py::arg("sf_method") =
            agnsf::esf::SFMethod::SecondOrder,
        py::arg("time") = "time",
        py::arg("value") = "value",
        py::arg("error") = "error",
        py::arg("uncertainty") =
            agnsf::esf::UncertaintyConfig(),
        py::arg("redshift") = 0.0,
        R"pbdoc(
Calculate the aggregated ensemble structure function from a list of
light-curve files (CSV or FITS).

method:    EnsembleMethod combination method
sf_method: SFMethod per-curve estimator
uncertainty: UncertaintyConfig (see ensemble_sf())
redshift: source redshift z applied to all light curves
        )pbdoc"
    );

    m.def(
        "ensemble_sf_from_path_list",
        [](
            const std::string& path_list_file,
            const agnsf::esf::LagBins& bins,
            agnsf::esf::SFEnsembleCalculator::Method method,
            agnsf::esf::SFMethod sf_method,
            const std::string& time,
            const std::string& value,
            const std::string& error,
            const agnsf::esf::UncertaintyConfig& uncertainty
        )
        {
            agnsf::io::ColumnNames columns;
            columns.time = time;
            columns.value = value;
            columns.error = error;

            return agnsf::esf::ensemble_sf_from_path_list(
                path_list_file,
                bins,
                sf_method,
                method,
                columns,
                uncertainty
            );
        },
        py::arg("path_list_file"),
        py::arg("bins"),
        py::arg("method") =
            agnsf::esf::SFEnsembleCalculator::Method::SqrtMeanSquared,
        py::arg("sf_method") =
            agnsf::esf::SFMethod::SecondOrder,
        py::arg("time") = "time",
        py::arg("value") = "value",
        py::arg("error") = "error",
        py::arg("uncertainty") =
            agnsf::esf::UncertaintyConfig(),
        R"pbdoc(
Calculate the aggregated ensemble structure function from a text file
that lists the light-curve paths (one per line).
        )pbdoc"
    );

    m.def(
        "write_sf_result",
        [](
            const std::string& path,
            const agnsf::esf::LagBins& bins,
            const agnsf::esf::SFResult& result
        )
        {
            agnsf::esf::write_sf_result(path, bins, result);
        },
        py::arg("path"),
        py::arg("bins"),
        py::arg("result"),
        R"pbdoc(
Write an SFResult to a simple text file.

Columns:

  # lag count sf_squared sf
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