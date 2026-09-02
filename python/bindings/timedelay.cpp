#include "common.hpp"
#include <optimization/optimization.hpp>
#include <timedelay/cross_correlation.hpp>
#include <timedelay/fr_rss.hpp>
#include <timedelay/transfer_function.hpp>

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


    // ------------------------------------------------------------------
    // Transfer-function fitting
    // ------------------------------------------------------------------
    py::enum_<agnsf::timedelay::TransferFunctionShape>(
        m,
        "TransferFunctionShape"
    )
    .value(
        "Gaussian",
        agnsf::timedelay::TransferFunctionShape::Gaussian
    )
    .value(
        "TopHat",
        agnsf::timedelay::TransferFunctionShape::TopHat
    );

    py::enum_<agnsf::optimization::Algorithm>(
        m,
        "OptimizationAlgorithm"
    )
    .value(
        "BoundedBobyqa",
        agnsf::optimization::Algorithm::BoundedBobyqa
    )
    .value(
        "NelderMead",
        agnsf::optimization::Algorithm::NelderMead
    );

    py::class_<agnsf::optimization::Parameter>(m, "FitParameter")
    .def(
        py::init<
            double,
            double,
            double,
            bool
        >(),
        py::arg("value") = 0.0,
        py::arg("lower") = 0.0,
        py::arg("upper") = 0.0,
        py::arg("fixed") = false
    )
    .def_readwrite(
        "value",
        &agnsf::optimization::Parameter::value
    )
    .def_readwrite(
        "lower",
        &agnsf::optimization::Parameter::lower
    )
    .def_readwrite(
        "upper",
        &agnsf::optimization::Parameter::upper
    )
    .def_readwrite(
        "fixed",
        &agnsf::optimization::Parameter::fixed
    );

    py::class_<agnsf::timedelay::TransferFunctionResult>(
        m,
        "TransferFunctionResult"
    )
    .def_readonly(
        "converged",
        &agnsf::timedelay::TransferFunctionResult::converged
    )
    .def_readonly(
        "offset",
        &agnsf::timedelay::TransferFunctionResult::offset
    )
    .def_readonly(
        "amplitude",
        &agnsf::timedelay::TransferFunctionResult::amplitude
    )
    .def_readonly(
        "lag",
        &agnsf::timedelay::TransferFunctionResult::lag
    )
    .def_readonly(
        "width",
        &agnsf::timedelay::TransferFunctionResult::width
    )
    .def_readonly(
        "chi2",
        &agnsf::timedelay::TransferFunctionResult::chi2
    )
    .def_readonly(
        "evaluations",
        &agnsf::timedelay::TransferFunctionResult::evaluations
    )
    .def_readonly(
        "n_valid_points",
        &agnsf::timedelay::TransferFunctionResult::n_valid_points
    )
    .def_readonly(
        "message",
        &agnsf::timedelay::TransferFunctionResult::message
    );

    m.def(
        "default_transfer_function_parameters",
        [](
            const Float64Array& time1,
            const Float64Array& value1,
            const Float64Array& error1,
            const Float64Array& time2,
            const Float64Array& value2,
            const Float64Array& error2,
            agnsf::timedelay::TransferFunctionShape shape,
            double grid_step
        )
        {
            const agnsf::LightCurve continuum =
                make_light_curve(time1, value1, error1);

            const agnsf::LightCurve response =
                make_light_curve(time2, value2, error2);

            agnsf::timedelay::TransferFunctionConfig config;
            config.shape = shape;
            config.grid_step = grid_step;

            return agnsf::timedelay::
                default_transfer_function_parameters(
                    continuum,
                    response,
                    config
                );
        },
        py::arg("time1"),
        py::arg("value1"),
        py::arg("error1"),
        py::arg("time2"),
        py::arg("value2"),
        py::arg("error2"),
        py::arg("shape") =
            agnsf::timedelay::TransferFunctionShape::Gaussian,
        py::arg("grid_step") = 0.0,
        R"pbdoc(
Data-driven initial values and box bounds for the four transfer
function parameters [offset, amplitude, lag, width].
        )pbdoc"
    );

    m.def(
        "fit_transfer_function",
        [](
            const Float64Array& time1,
            const Float64Array& value1,
            const Float64Array& error1,
            const Float64Array& time2,
            const Float64Array& value2,
            const Float64Array& error2,
            const std::vector<agnsf::optimization::Parameter>& parameters,
            agnsf::timedelay::TransferFunctionShape shape,
            double grid_step,
            std::size_t lag_restarts,
            agnsf::optimization::Algorithm algorithm,
            double xtol_rel,
            double ftol_rel,
            std::size_t max_evaluations
        )
        {
            const agnsf::LightCurve continuum =
                make_light_curve(time1, value1, error1);

            const agnsf::LightCurve response =
                make_light_curve(time2, value2, error2);

            agnsf::timedelay::TransferFunctionConfig config;
            config.shape = shape;
            config.grid_step = grid_step;
            config.lag_restarts = lag_restarts;

            agnsf::optimization::Options options;
            options.algorithm = algorithm;
            options.xtol_rel = xtol_rel;
            options.ftol_rel = ftol_rel;
            options.max_evaluations = max_evaluations;

            return agnsf::timedelay::fit_transfer_function(
                continuum,
                response,
                parameters,
                config,
                options
            );
        },
        py::arg("time1"),
        py::arg("value1"),
        py::arg("error1"),
        py::arg("time2"),
        py::arg("value2"),
        py::arg("error2"),
        py::arg("parameters"),
        py::arg("shape") =
            agnsf::timedelay::TransferFunctionShape::Gaussian,
        py::arg("grid_step") = 0.0,
        py::arg("lag_restarts") = std::size_t{0},
        py::arg("algorithm") =
            agnsf::optimization::Algorithm::BoundedBobyqa,
        py::arg("xtol_rel") = 1e-6,
        py::arg("ftol_rel") = 1e-8,
        py::arg("max_evaluations") = std::size_t{10000},
        R"pbdoc(
Fit a parametric transfer function Psi(tau) by chi^2 minimization.

Model: R(t) = offset + amplitude * (C * Psi)(t), where Psi is a
normalized Gaussian or top-hat centered at `lag` with `width`.

parameters: four FitParameter entries [offset, amplitude, lag,
            width] giving initial values and box bounds; lower ==
            upper (or fixed=True) pins a parameter.

The linear parameters (offset, amplitude) are profiled
analytically; only (lag, width) are searched. Returns a
TransferFunctionResult with the best parameters and chi^2.
        )pbdoc"
    );

    m.def(
        "transfer_function_curve",
        [](
            const std::vector<double>& taus,
            agnsf::timedelay::TransferFunctionShape shape,
            double lag,
            double width
        )
        {
            return agnsf::timedelay::transfer_function_curve(
                taus,
                shape,
                lag,
                width
            );
        },
        py::arg("taus"),
        py::arg("shape"),
        py::arg("lag"),
        py::arg("width"),
        R"pbdoc(
Evaluate the normalized transfer function Psi(tau) on a lag grid.
        )pbdoc"
    );

    m.def(
        "transfer_function_model_response",
        [](
            const Float64Array& time1,
            const Float64Array& value1,
            const Float64Array& error1,
            const std::vector<double>& response_times,
            agnsf::timedelay::TransferFunctionShape shape,
            double offset,
            double amplitude,
            double lag,
            double width,
            double grid_step
        )
        {
            const agnsf::LightCurve continuum =
                make_light_curve(time1, value1, error1);

            return agnsf::timedelay::
                transfer_function_model_response(
                    continuum,
                    response_times,
                    shape,
                    offset,
                    amplitude,
                    lag,
                    width,
                    grid_step
                );
        },
        py::arg("time1"),
        py::arg("value1"),
        py::arg("error1"),
        py::arg("response_times"),
        py::arg("shape"),
        py::arg("offset"),
        py::arg("amplitude"),
        py::arg("lag"),
        py::arg("width"),
        py::arg("grid_step") = 0.0,
        R"pbdoc(
Model response R(t) = offset + amplitude * (C * Psi)(t) at the
given response times for fixed transfer-function parameters.
        )pbdoc"
    );
}
