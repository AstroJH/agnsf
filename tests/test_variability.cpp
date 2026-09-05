#include <cassert>
#include <cmath>
#include <vector>

#include <core/light_curve.hpp>
#include <variability/variability.hpp>

namespace {

bool close(
    double a,
    double b,
    double rtol = 1e-6,
    double atol = 1e-8
)
{
    return std::abs(a - b) <= atol + rtol * std::abs(b);
}


agnsf::LightCurve make_curve(
    const std::vector<double>& values,
    const std::vector<double>& errors
)
{
    std::vector<double> time;

    for (std::size_t i = 0; i < values.size(); ++i) {
        time.push_back(static_cast<double>(i));
    }

    return agnsf::LightCurve(time, values, errors);
}


void test_constant_curve()
{
    // Constant source: all variability must vanish, only noise.
    const agnsf::LightCurve curve =
        make_curve({10.0, 10.0, 10.0, 10.0}, {1.0, 1.0, 1.0, 1.0});

    const agnsf::variability::Statistics s =
        agnsf::variability::measure(curve.view());

    assert(s.valid);
    assert(s.n == 4);
    assert(close(s.mean, 10.0));
    assert(close(s.stddev, 0.0));
    assert(close(s.peak_to_peak, 0.0));
    assert(close(s.sigma_m, 0.0));
    assert(close(s.fvar, 0.0));
    assert(close(s.nxs, -0.01)); // (S2 - <e^2>)/mean^2 = -1/100
    assert(close(s.xs, -1.0));   // S2 - <e^2> = 0 - 1
    assert(close(s.chi2, 0.0));
    assert(close(s.chi2_q, 1.0)); // Q(a, 0) = 1
    assert(!s.fvar_uncertainty.estimated());

    // Weighted mean is only computed on request (Options::weighted).
    assert(std::isnan(s.weighted_mean));

    agnsf::variability::Options weighted;
    weighted.weighted = true;

    const agnsf::variability::Statistics ws =
        agnsf::variability::measure(curve.view(), weighted);

    assert(close(ws.weighted_mean, 10.0));
    assert(close(ws.weighted_mean_error, 0.5));
}


void test_sample_statistics()
{
    // Linear ramp 0..9: analytic mean / variance / std / PP.
    const agnsf::LightCurve curve =
        make_curve({0, 1, 2, 3, 4, 5, 6, 7, 8, 9}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0});

    const agnsf::variability::Statistics s =
        agnsf::variability::measure(curve.view());

    assert(close(s.mean, 4.5));
    assert(close(s.stddev, std::sqrt(82.5 / 9.0)));
    assert(close(s.peak_to_peak, 9.0));
    assert(close(s.peak_to_peak_noise_corrected, 9.0)); // no noise

    // Von Neumann ratio: sum of successive-diff^2 / sum(x - mean)^2.
    assert(close(s.von_neumann, 9.0 / 82.5));

    // No positive errors: chi^2 / weighted mean are undefined.
    assert(std::isnan(s.chi2));
    assert(std::isnan(s.weighted_mean));

    // Constant curve: von Neumann ratio undefined (zero denominator).
    const agnsf::LightCurve flat =
        make_curve({3.0, 3.0, 3.0}, {0.0, 0.0, 0.0});
    const agnsf::variability::Statistics fs =
        agnsf::variability::measure(flat.view());
    assert(std::isnan(fs.von_neumann));
}


void test_weighted_mean()
{
    const agnsf::LightCurve curve =
        make_curve({0.0, 10.0}, {1.0, 3.0});

    agnsf::variability::Options options;
    options.weighted = true;

    const agnsf::variability::Statistics s =
        agnsf::variability::measure(curve.view(), options);

    // weights 1 and 1/9 -> weighted mean = (10/9)/(10/9) = 1
    assert(close(s.weighted_mean, 1.0));
    // error = 1 / sqrt(10/9) = 3 / sqrt(10)
    assert(close(s.weighted_mean_error, 3.0 / std::sqrt(10.0)));
}


void test_fvar_and_excess_variance()
{
    // values [2, 4, 6], errors 1: xbar = 4, S2 = 4, <e^2> = 1
    const agnsf::LightCurve curve =
        make_curve({2.0, 4.0, 6.0}, {1.0, 1.0, 1.0});

    const agnsf::variability::Statistics s =
        agnsf::variability::measure(curve.view());

    assert(close(s.xs, 3.0));
    assert(close(s.nxs, 3.0 / 16.0));
    assert(close(s.fvar, std::sqrt(3.0 / 16.0)));

    // Vaughan et al. (2003) analytic error.
    const double eps2 = 1.0;
    const double xbar2 = 16.0;
    const double nn = 3.0;
    const double fvar = s.fvar;

    const double expected_error = std::sqrt(
        eps2 / (nn * xbar2) +
        0.5 / nn * std::pow(eps2 / (xbar2 * fvar), 2.0)
    );

    assert(s.fvar_uncertainty.estimated());
    assert(close(s.fvar_uncertainty.lower, fvar - expected_error, 1e-4, 1e-8));
    assert(close(s.fvar_uncertainty.upper, fvar + expected_error, 1e-4, 1e-8));

    // chi^2 = (4 + 0 + 4)/1 = 8, dof = 2, Q = Q(1, 4) = exp(-4).
    assert(close(s.chi2, 8.0));
    assert(close(s.chi2_dof, 2.0));
    assert(close(s.chi2_q, std::exp(-4.0), 1e-6, 1e-10));
}


void test_chi2_survival_identity()
{
    // Two points -> dof = 1 -> Q(0.5, chi2/2) = erfc(sqrt(chi2/2)).
    const agnsf::LightCurve curve =
        make_curve({0.0, 1.0}, {1.0, 1.0});

    const agnsf::variability::Statistics s =
        agnsf::variability::measure(curve.view());

    // mean = 0.5 -> chi2 = 0.25 + 0.25 = 0.5
    assert(close(s.chi2, 0.5));
    assert(close(s.chi2_q, std::erfc(std::sqrt(0.25)), 1e-9, 1e-12));
}


void test_noise_corrected_peak_to_peak()
{
    const agnsf::LightCurve curve =
        make_curve({0.0, 3.0}, {1.0, 1.0});

    const agnsf::variability::Statistics s =
        agnsf::variability::measure(curve.view());

    // range = 3, 2*<e^2> = 2 -> sqrt(9 - 2) = sqrt(7)
    assert(close(s.peak_to_peak, 3.0));
    assert(close(s.peak_to_peak_noise_corrected, std::sqrt(7.0)));
}


void test_err_sys()
{
    // values [0, 1, 2] with no measurement errors but a systematic
    // floor of 1: sample variance 1 is fully absorbed -> sigma_m = 0.
    const agnsf::LightCurve curve =
        make_curve({0.0, 1.0, 2.0}, {0.0, 0.0, 0.0});

    agnsf::variability::Options options;
    options.err_sys = 1.0;

    const agnsf::variability::Statistics s =
        agnsf::variability::measure(curve.view(), options);

    assert(close(s.stddev, 1.0));
    assert(close(s.sigma_m, 0.0));
    assert(close(s.nxs, 0.0)); // S2 - e^2 = 0
}


void test_too_few_points()
{
    const agnsf::LightCurve curve =
        make_curve({1.0}, {0.1});

    const agnsf::variability::Statistics s =
        agnsf::variability::measure(curve.view());

    assert(!s.valid);
}

} // namespace


int main()
{
    test_constant_curve();
    test_sample_statistics();
    test_weighted_mean();
    test_fvar_and_excess_variance();
    test_chi2_survival_identity();
    test_noise_corrected_peak_to_peak();
    test_err_sys();
    test_too_few_points();

    return 0;
}
