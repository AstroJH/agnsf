#include "common.hpp"
#include <variability/variability.hpp>

using namespace agnsf::python;


// ------------------------------------------------------------------
// Single-curve variability statistics
// ------------------------------------------------------------------

void bind_variability(py::module_& m)
{
    py::class_<agnsf::variability::Statistics>(
        m,
        "VariabilityResult"
    )
    .def_readonly(
        "n",
        &agnsf::variability::Statistics::n
    )
    .def_readonly(
        "valid",
        &agnsf::variability::Statistics::valid
    )
    .def_readonly(
        "mean",
        &agnsf::variability::Statistics::mean
    )
    .def_readonly(
        "weighted_mean",
        &agnsf::variability::Statistics::weighted_mean
    )
    .def_readonly(
        "weighted_mean_error",
        &agnsf::variability::Statistics::weighted_mean_error
    )
    .def_readonly(
        "stddev",
        &agnsf::variability::Statistics::stddev
    )
    .def_readonly(
        "stddev_error",
        &agnsf::variability::Statistics::stddev_error
    )
    .def_readonly(
        "peak_to_peak",
        &agnsf::variability::Statistics::peak_to_peak
    )
    .def_readonly(
        "peak_to_peak_noise_corrected",
        &agnsf::variability::Statistics::peak_to_peak_noise_corrected
    )
    .def_readonly(
        "sigma_m",
        &agnsf::variability::Statistics::sigma_m
    )
    .def_readonly(
        "fvar",
        &agnsf::variability::Statistics::fvar
    )
    .def_readonly(
        "fvar_uncertainty",
        &agnsf::variability::Statistics::fvar_uncertainty
    )
    .def_readonly(
        "nxs",
        &agnsf::variability::Statistics::nxs
    )
    .def_readonly(
        "nxs_uncertainty",
        &agnsf::variability::Statistics::nxs_uncertainty
    )
    .def_readonly(
        "xs",
        &agnsf::variability::Statistics::xs
    )
    .def_readonly(
        "xs_uncertainty",
        &agnsf::variability::Statistics::xs_uncertainty
    )
    .def_readonly(
        "chi2",
        &agnsf::variability::Statistics::chi2
    )
    .def_readonly(
        "chi2_dof",
        &agnsf::variability::Statistics::chi2_dof
    )
    .def_readonly(
        "chi2_q",
        &agnsf::variability::Statistics::chi2_q
    )
    .def_readonly(
        "von_neumann",
        &agnsf::variability::Statistics::von_neumann
    );

    m.def(
        "variability_measure",
        [](
            const Float64Array& time,
            const Float64Array& value,
            const Float64Array& error,
            double err_sys,
            bool weighted
        )
        {
            const agnsf::LightCurveView curve =
                make_light_curve_view(time, value, error);

            agnsf::variability::Options options;
            options.err_sys = err_sys;
            options.weighted = weighted;

            return agnsf::variability::measure(curve, options);
        },
        py::arg("time"),
        py::arg("value"),
        py::arg("error"),
        py::arg("err_sys") = 0.0,
        py::arg("weighted") = false,
        R"pbdoc(
Compute single-curve variability statistics.

Returns a VariabilityResult with the weighted/plain mean, sample
standard deviation, peak-to-peak amplitude, intrinsic sigma_m,
fractional variability amplitude F_var, normalized and
unnormalized excess variance (with Vaughan et al. 2003 analytic
uncertainties), the chi^2 variability significance Q, and the von
Neumann ratio.

err_sys: systematic error floor added in quadrature when the noise
            contribution is subtracted.
weighted: use the inverse-variance weighted mean inside the
            variance estimates.
        )pbdoc"
    );
}
