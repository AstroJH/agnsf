#include <cassert>
#include <cmath>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#include <core/light_curve.hpp>
#include <timedelay/transfer_function.hpp>

namespace {

constexpr double kPi = 3.14159265358979323846;

bool close(
    double a,
    double b,
    double rtol = 1e-4,
    double atol = 1e-4
)
{
    return std::abs(a - b) <= atol + rtol * std::abs(b);
}


agnsf::LightCurve make_continuum()
{
    // Smooth analytic continuum; the reference convolution evaluates
    // it directly, so no interpolation is involved on the test side.
    std::vector<double> t;
    std::vector<double> v;
    std::vector<double> e;

    for (int i = 0; i <= 200; ++i) {
        const double x = static_cast<double>(i);
        t.push_back(x);
        v.push_back(
            2.0
            + std::sin(2.0 * kPi * x / 50.0)
            + 0.5 * std::sin(2.0 * kPi * x / 23.0)
            + 0.3 * std::cos(2.0 * kPi * x / 11.0)
        );
        e.push_back(0.01);
    }

    return agnsf::LightCurve(t, v, e);
}


double continuum_analytic(double s)
{
    return
        2.0
        + std::sin(2.0 * kPi * s / 50.0)
        + 0.5 * std::sin(2.0 * kPi * s / 23.0)
        + 0.3 * std::cos(2.0 * kPi * s / 11.0);
}


double psi_reference(
    double tau,
    agnsf::timedelay::TransferFunctionShape shape,
    double lag,
    double width
)
{
    const double z = tau - lag;

    if (shape == agnsf::timedelay::TransferFunctionShape::TopHat) {
        return std::abs(z) <= width ? 1.0 / (2.0 * width) : 0.0;
    }

    const double u = z / width;
    return std::exp(-0.5 * u * u) /
        (std::sqrt(2.0 * kPi) * width);
}


// Independent reference response: dense-grid convolution using the
// analytic continuum (step 0.02, rectangle rule).
double reference_response(
    double t,
    agnsf::timedelay::TransferFunctionShape shape,
    double offset,
    double amplitude,
    double lag,
    double width
)
{
    const double half =
        shape == agnsf::timedelay::TransferFunctionShape::TopHat
            ? width
            : 4.0 * width;

    constexpr double step = 0.02;

    double integral = 0.0;

    for (double s = t - lag - half;
         s <= t - lag + half + 1e-9;
         s += step) {

        const double tau = t - s;
        integral +=
            continuum_analytic(s) *
            psi_reference(tau, shape, lag, width) *
            step;
    }

    return offset + amplitude * integral;
}


agnsf::LightCurve make_response(
    agnsf::timedelay::TransferFunctionShape shape,
    double offset,
    double amplitude,
    double lag,
    double width,
    double noise,
    std::vector<double>* observed_times = nullptr
)
{
    const double half =
        shape == agnsf::timedelay::TransferFunctionShape::TopHat
            ? width
            : 4.0 * width;

    // Keep only times whose convolution window lies well inside the
    // continuum coverage [0, 200].
    const int margin = static_cast<int>(std::ceil(half + lag)) + 5;
    const int t_start = margin;
    const int t_end = 200 - margin;

    std::mt19937 rng(42);
    std::normal_distribution<double> normal(0.0, noise);

    std::vector<double> t;
    std::vector<double> v;
    std::vector<double> e;

    for (int i = t_start; i <= t_end; ++i) {
        const double x = static_cast<double>(i);
        t.push_back(x);

        const double model = reference_response(
            x, shape, offset, amplitude, lag, width
        );

        v.push_back(model + normal(rng));
        e.push_back(noise);

        if (observed_times != nullptr) {
            observed_times->push_back(x);
        }
    }

    return agnsf::LightCurve(t, v, e);
}


void test_gaussian_recovery()
{
    const double true_offset = 0.4;
    const double true_amplitude = 1.5;
    const double true_lag = 12.0;
    const double true_width = 3.0;

    const agnsf::LightCurve continuum = make_continuum();

    std::vector<double> times;

    const agnsf::LightCurve response = make_response(
        agnsf::timedelay::TransferFunctionShape::Gaussian,
        true_offset,
        true_amplitude,
        true_lag,
        true_width,
        0.01,
        &times
    );

    agnsf::timedelay::TransferFunctionConfig config;
    config.shape = agnsf::timedelay::TransferFunctionShape::Gaussian;
    config.grid_step = 0.1;

    std::vector<agnsf::optimization::Parameter> parameters = {
        {0.0, -5.0, 5.0},    // offset
        {1.0, 0.0, 10.0},    // amplitude
        {5.0, 0.0, 40.0},    // lag
        {2.0, 0.2, 10.0},    // width
    };

    const agnsf::timedelay::TransferFunctionResult result =
        agnsf::timedelay::fit_transfer_function(
            continuum,
            response,
            parameters,
            config
        );

    assert(result.converged);
    assert(result.n_valid_points == response.size());

    assert(close(result.lag, true_lag, 0.05, 0.3));
    assert(close(result.width, true_width, 0.05, 0.4));
    assert(close(result.amplitude, true_amplitude, 0.1, 0.1));
    assert(close(result.offset, true_offset, 0.1, 0.1));

    // End-to-end: the fitted model response must reproduce the
    // injected noiseless response.
    const std::vector<double> model =
        agnsf::timedelay::transfer_function_model_response(
            continuum,
            times,
            agnsf::timedelay::TransferFunctionShape::Gaussian,
            result.offset,
            result.amplitude,
            result.lag,
            result.width,
            config.grid_step
        );

    for (std::size_t i = 0; i < times.size(); ++i) {
        const double expected = reference_response(
            times[i],
            agnsf::timedelay::TransferFunctionShape::Gaussian,
            true_offset,
            true_amplitude,
            true_lag,
            true_width
        );

        assert(std::abs(model[i] - expected) < 0.05);
    }
}


void test_tophat_recovery()
{
    const double true_offset = 0.2;
    const double true_amplitude = 2.0;
    const double true_lag = 10.0;
    const double true_width = 2.0;

    const agnsf::LightCurve continuum = make_continuum();

    const agnsf::LightCurve response = make_response(
        agnsf::timedelay::TransferFunctionShape::TopHat,
        true_offset,
        true_amplitude,
        true_lag,
        true_width,
        0.01
    );

    agnsf::timedelay::TransferFunctionConfig config;
    config.shape = agnsf::timedelay::TransferFunctionShape::TopHat;
    config.grid_step = 0.1;

    std::vector<agnsf::optimization::Parameter> parameters = {
        {0.0, -5.0, 5.0},
        {1.0, 0.0, 10.0},
        {0.0, 0.0, 40.0},
        {1.0, 0.2, 10.0},
    };

    const agnsf::timedelay::TransferFunctionResult result =
        agnsf::timedelay::fit_transfer_function(
            continuum,
            response,
            parameters,
            config
        );

    assert(result.converged);
    assert(close(result.lag, true_lag, 0.05, 0.4));
    assert(close(result.width, true_width, 0.05, 0.5));
    assert(close(result.amplitude, true_amplitude, 0.1, 0.15));
    assert(close(result.offset, true_offset, 0.1, 0.15));
}


void test_curve_normalization()
{
    const double lag = 10.0;
    const double width = 2.0;

    const std::vector<double> taus = [] {
        std::vector<double> out;
        for (double tau = -20.0; tau <= 40.0 + 1e-9; tau += 0.01) {
            out.push_back(tau);
        }
        return out;
    }();

    const std::vector<double> gaussian =
        agnsf::timedelay::transfer_function_curve(
            taus,
            agnsf::timedelay::TransferFunctionShape::Gaussian,
            lag,
            width
        );

    const std::vector<double> tophat =
        agnsf::timedelay::transfer_function_curve(
            taus,
            agnsf::timedelay::TransferFunctionShape::TopHat,
            lag,
            width
        );

    auto integral = [&taus](const std::vector<double>& values) {
        double sum = 0.0;

        for (std::size_t i = 1; i < values.size(); ++i) {
            sum += 0.5 * (values[i - 1] + values[i]) *
                (taus[i] - taus[i - 1]);
        }

        return sum;
    };

    assert(close(integral(gaussian), 1.0, 1e-3, 1e-3));
    assert(close(integral(tophat), 1.0, 1e-3, 1e-3));

    // Evaluate at a single lag for the peak value.
    const double peak = agnsf::timedelay::evaluate_transfer_function(
        lag,
        agnsf::timedelay::TransferFunctionShape::Gaussian,
        lag,
        width
    );
    assert(close(peak, 1.0 / (width * std::sqrt(2.0 * kPi)), 1e-6, 1e-6));
}


void test_fixed_nonlinear_parameters()
{
    const double true_offset = 0.3;
    const double true_amplitude = 1.2;
    const double true_lag = 8.0;
    const double true_width = 2.5;

    const agnsf::LightCurve continuum = make_continuum();

    const agnsf::LightCurve response = make_response(
        agnsf::timedelay::TransferFunctionShape::Gaussian,
        true_offset,
        true_amplitude,
        true_lag,
        true_width,
        0.005
    );

    agnsf::timedelay::TransferFunctionConfig config;
    config.shape = agnsf::timedelay::TransferFunctionShape::Gaussian;
    config.grid_step = 0.1;

    // Pin lag and width at their true values; recover offset and
    // amplitude only (linear parameters).
    std::vector<agnsf::optimization::Parameter> parameters = {
        {0.0, -5.0, 5.0},
        {0.5, 0.0, 10.0},
        {true_lag, true_lag, true_lag},
        {true_width, true_width, true_width},
    };

    const agnsf::timedelay::TransferFunctionResult result =
        agnsf::timedelay::fit_transfer_function(
            continuum,
            response,
            parameters,
            config
        );

    assert(result.converged);
    assert(close(result.offset, true_offset, 0.1, 0.05));
    assert(close(result.amplitude, true_amplitude, 0.1, 0.05));
    assert(result.lag == true_lag);
    assert(result.width == true_width);
}


void test_default_parameters()
{
    const agnsf::LightCurve continuum = make_continuum();

    std::vector<double> times;

    const agnsf::LightCurve response = make_response(
        agnsf::timedelay::TransferFunctionShape::Gaussian,
        0.4,
        1.5,
        12.0,
        3.0,
        0.01,
        &times
    );

    const std::vector<agnsf::optimization::Parameter> parameters =
        agnsf::timedelay::default_transfer_function_parameters(
            continuum,
            response
        );

    assert(parameters.size() == 4);
    assert(parameters[3].lower > 0.0);
    assert(parameters[1].lower == 0.0);
    assert(parameters[2].lower <= 12.0);
    assert(parameters[2].upper >= 12.0);

    // The defaults must be a valid input for the fit (auto lag scan).
    agnsf::timedelay::TransferFunctionConfig config;

    const agnsf::timedelay::TransferFunctionResult result =
        agnsf::timedelay::fit_transfer_function(
            continuum,
            response,
            parameters,
            config
        );

    assert(result.converged);
    assert(close(result.lag, 12.0, 0.1, 2.0));
}


void test_exceptions()
{
    const agnsf::LightCurve continuum = make_continuum();

    const agnsf::LightCurve response = make_response(
        agnsf::timedelay::TransferFunctionShape::Gaussian,
        0.4,
        1.5,
        12.0,
        3.0,
        0.01
    );

    agnsf::timedelay::TransferFunctionConfig config;

    const std::vector<agnsf::optimization::Parameter> valid = {
        {0.0, -5.0, 5.0},
        {1.0, 0.0, 10.0},
        {5.0, 0.0, 40.0},
        {2.0, 0.2, 10.0},
    };

    bool threw = false;

    // Wrong parameter count.
    try {
        agnsf::timedelay::fit_transfer_function(
            continuum, response,
            {valid[0], valid[1], valid[2]},
            config
        );
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);

    // width lower bound not positive.
    threw = false;

    try {
        std::vector<agnsf::optimization::Parameter> bad = valid;
        bad[3].lower = 0.0;

        agnsf::timedelay::fit_transfer_function(
            continuum, response, bad, config
        );
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);

    // lower > upper.
    threw = false;

    try {
        std::vector<agnsf::optimization::Parameter> bad = valid;
        bad[2].lower = 40.0;
        bad[2].upper = 0.0;

        agnsf::timedelay::fit_transfer_function(
            continuum, response, bad, config
        );
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);

    // Negative grid step.
    threw = false;

    try {
        agnsf::timedelay::TransferFunctionConfig bad_config;
        bad_config.grid_step = -1.0;

        agnsf::timedelay::fit_transfer_function(
            continuum, response, valid, bad_config
        );
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);

    // Degenerate continuum.
    threw = false;

    try {
        const agnsf::LightCurve flat(
            {1.0, 1.0},
            {0.0, 1.0},
            {0.1, 0.1}
        );

        agnsf::timedelay::fit_transfer_function(
            flat, response, valid, config
        );
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);

    // Non-finite data.
    threw = false;

    try {
        const agnsf::LightCurve bad(
            {0.0, 1.0, 2.0},
            {0.0, std::numeric_limits<double>::quiet_NaN(), 2.0},
            {0.1, 0.1, 0.1}
        );

        agnsf::timedelay::fit_transfer_function(
            bad, response, valid, config
        );
    } catch (const std::invalid_argument&) {
        threw = true;
    }
    assert(threw);
}

} // namespace


int main()
{
    test_gaussian_recovery();
    test_tophat_recovery();
    test_curve_normalization();
    test_fixed_nonlinear_parameters();
    test_default_parameters();
    test_exceptions();

    return 0;
}
