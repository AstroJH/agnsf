#include "common.hpp"
#include <esf/pooled_esf_calculator.hpp>
#include <esf/sf_ensemble_calculator.hpp>
#include <esf/sf_io.hpp>
#include <esf/sf_result.hpp>

void bind_esf(py::module_& m)
{
    m.def(
            "pooled_sf",
            [](
                const py::list& times,
                const py::list& values,
                const py::list& errors,
                const agnsf::esf::LagBins& bins,
                agnsf::esf::SFMethod method,
                const agnsf::UncertaintyConfig& uncertainty,
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
                agnsf::UncertaintyConfig(),
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
                const agnsf::UncertaintyConfig& uncertainty,
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
                agnsf::UncertaintyConfig(),
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
            "sf_from_file",
            [](
                const std::string& path,
                const agnsf::esf::LagBins& bins,
                agnsf::esf::SFMethod method,
                const std::string& time,
                const std::string& value,
                const std::string& error,
                const agnsf::UncertaintyConfig& uncertainty,
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
                agnsf::UncertaintyConfig(),
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
                const agnsf::UncertaintyConfig& uncertainty,
                const py::object& redshift
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
    
                if (py::isinstance<py::float_>(redshift) ||
                    py::isinstance<py::int_>(redshift)) {
    
                    // Scalar: the same z for every light curve.
                    return agnsf::esf::pooled_sf_from_files(
                        path_vector,
                        bins,
                        method,
                        columns,
                        uncertainty,
                        redshift.cast<double>()
                    );
                }
    
                // Sequence: one redshift per light curve.
                const std::vector<double> redshifts =
                    redshift.cast<std::vector<double>>();
    
                return agnsf::esf::pooled_sf_from_files(
                    path_vector,
                    redshifts,
                    bins,
                    method,
                    columns,
                    uncertainty
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
                agnsf::UncertaintyConfig(),
            py::arg("redshift") = 0.0,
            R"pbdoc(
    Calculate the pooled ensemble structure function from a list of
    light-curve files (CSV or FITS).
    
    redshift: source redshift z — a scalar applies the same z to all
                light curves; a sequence provides one z per light curve.
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
                const agnsf::UncertaintyConfig& uncertainty
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
                agnsf::UncertaintyConfig(),
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
                const agnsf::UncertaintyConfig& uncertainty,
                const py::object& redshift
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
    
                if (py::isinstance<py::float_>(redshift) ||
                    py::isinstance<py::int_>(redshift)) {
    
                    // Scalar: the same z for every light curve.
                    return agnsf::esf::ensemble_sf_from_files(
                        path_vector,
                        bins,
                        sf_method,
                        method,
                        columns,
                        uncertainty,
                        redshift.cast<double>()
                    );
                }
    
                // Sequence: one redshift per light curve.
                const std::vector<double> redshifts =
                    redshift.cast<std::vector<double>>();
    
                return agnsf::esf::ensemble_sf_from_files(
                    path_vector,
                    redshifts,
                    bins,
                    sf_method,
                    method,
                    columns,
                    uncertainty
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
                agnsf::UncertaintyConfig(),
            py::arg("redshift") = 0.0,
            R"pbdoc(
    Calculate the aggregated ensemble structure function from a list of
    light-curve files (CSV or FITS).
    
    method:    EnsembleMethod combination method
    sf_method: SFMethod per-curve estimator
    uncertainty: UncertaintyConfig (see ensemble_sf())
    redshift: source redshift z — a scalar applies the same z to all
                light curves; a sequence provides one z per light curve.
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
                const agnsf::UncertaintyConfig& uncertainty
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
                agnsf::UncertaintyConfig(),
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
        // Time-delay (cross-correlation lag analysis)
        // ------------------------------------------------------------------
}
