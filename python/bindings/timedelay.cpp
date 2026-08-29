#include "common.hpp"
#include <timedelay/cross_correlation.hpp>
#include <timedelay/fr_rss.hpp>

using namespace agnsf::python;


// ------------------------------------------------------------------
// Time-delay (cross-correlation lag analysis)
// ------------------------------------------------------------------

void bind_timedelay(py::module_& m)
{
    py::enum_<agnsf::timedelay::CrossCorrelationMethod>(
            m,
            "CrossCorrelationMethod"
        )
            .value(
                "Dcf",
                agnsf::timedelay::CrossCorrelationMethod::Dcf
            )
            .value(
                "Iccf",
                agnsf::timedelay::CrossCorrelationMethod::Iccf
            );
    py::enum_<agnsf::timedelay::LagEstimate>(
            m,
            "LagEstimate"
        )
            .value(
                "Peak",
                agnsf::timedelay::LagEstimate::Peak
            )
            .value(
                "Centroid",
                agnsf::timedelay::LagEstimate::Centroid
            );
    py::class_<agnsf::timedelay::LagResult>(m, "LagResult")
            .def_readonly(
                "lag_peak",
                &agnsf::timedelay::LagResult::lag_peak
            )
            .def_readonly(
                "lag_centroid",
                &agnsf::timedelay::LagResult::lag_centroid
            )
            .def_readonly(
                "peak_value",
                &agnsf::timedelay::LagResult::peak_value
            )
            .def_readonly(
                "tau",
                &agnsf::timedelay::LagResult::tau
            )
            .def_readonly(
                "ccf",
                &agnsf::timedelay::LagResult::ccf
            )
            .def_readonly(
                "count",
                &agnsf::timedelay::LagResult::count
            );
    m.def(
            "cross_correlate",
            [](
                const Float64Array& time1,
                const Float64Array& value1,
                const Float64Array& error1,
                const Float64Array& time2,
                const Float64Array& value2,
                const Float64Array& error2,
                double grid_min,
                double grid_max,
                double grid_step,
                agnsf::timedelay::CrossCorrelationMethod method,
                double dcf_bin_width,
                double centroid_threshold,
                std::size_t min_overlap
            )
            {
                const agnsf::LightCurve continuum =
                    make_light_curve(time1, value1, error1);
    
                const agnsf::LightCurve response =
                    make_light_curve(time2, value2, error2);
    
                agnsf::timedelay::LagGrid grid;
                grid.min = grid_min;
                grid.max = grid_max;
                grid.step = grid_step;
    
                agnsf::timedelay::CrossCorrelationConfig config;
                config.method = method;
                config.dcf_bin_width = dcf_bin_width;
                config.centroid_threshold = centroid_threshold;
                config.min_overlap = min_overlap;
    
                return agnsf::timedelay::cross_correlate(
                    continuum,
                    response,
                    grid,
                    config
                );
            },
            py::arg("time1"),
            py::arg("value1"),
            py::arg("error1"),
            py::arg("time2"),
            py::arg("value2"),
            py::arg("error2"),
            py::arg("grid_min"),
            py::arg("grid_max"),
            py::arg("grid_step"),
            py::arg("method") =
                agnsf::timedelay::CrossCorrelationMethod::Dcf,
            py::arg("dcf_bin_width") = 1.0,
            py::arg("centroid_threshold") = 0.8,
            py::arg("min_overlap") = std::size_t{3},
            R"pbdoc(
    Cross-correlate a continuum light curve against a response light curve
    over a grid of trial lags.
    
    Positive lag means the response lags the continuum.
            )pbdoc"
        );
    m.def(
            "lag_uncertainty",
            [](
                const Float64Array& time1,
                const Float64Array& value1,
                const Float64Array& error1,
                const Float64Array& time2,
                const Float64Array& value2,
                const Float64Array& error2,
                double grid_min,
                double grid_max,
                double grid_step,
                agnsf::timedelay::LagEstimate estimate,
                agnsf::timedelay::CrossCorrelationMethod method,
                double dcf_bin_width,
                double centroid_threshold,
                std::size_t min_overlap,
                std::size_t n_realizations,
                std::uint32_t seed,
                bool flux_randomization,
                bool random_subset
            )
            {
                const agnsf::LightCurve continuum =
                    make_light_curve(time1, value1, error1);
    
                const agnsf::LightCurve response =
                    make_light_curve(time2, value2, error2);
    
                agnsf::timedelay::LagGrid grid;
                grid.min = grid_min;
                grid.max = grid_max;
                grid.step = grid_step;
    
                agnsf::timedelay::CrossCorrelationConfig config;
                config.method = method;
                config.dcf_bin_width = dcf_bin_width;
                config.centroid_threshold = centroid_threshold;
                config.min_overlap = min_overlap;
    
                agnsf::timedelay::FRRSSConfig fr_rss;
                fr_rss.n_realizations = n_realizations;
                fr_rss.seed = seed;
                fr_rss.flux_randomization = flux_randomization;
                fr_rss.random_subset = random_subset;
    
                return agnsf::timedelay::lag_uncertainty(
                    continuum,
                    response,
                    grid,
                    estimate,
                    config,
                    fr_rss
                );
            },
            py::arg("time1"),
            py::arg("value1"),
            py::arg("error1"),
            py::arg("time2"),
            py::arg("value2"),
            py::arg("error2"),
            py::arg("grid_min"),
            py::arg("grid_max"),
            py::arg("grid_step"),
            py::arg("estimate") =
                agnsf::timedelay::LagEstimate::Peak,
            py::arg("method") =
                agnsf::timedelay::CrossCorrelationMethod::Dcf,
            py::arg("dcf_bin_width") = 1.0,
            py::arg("centroid_threshold") = 0.8,
            py::arg("min_overlap") = std::size_t{3},
            py::arg("n_realizations") = std::size_t{1000},
            py::arg("seed") = std::uint32_t{0},
            py::arg("flux_randomization") = true,
            py::arg("random_subset") = true,
            R"pbdoc(
    FR/RSS Monte Carlo uncertainty of a lag estimate (Peterson et al. 1998).
            )pbdoc"
        );
}
