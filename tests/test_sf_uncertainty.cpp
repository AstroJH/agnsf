#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <core/light_curve.hpp>
#include <esf/sf_calculator.hpp>
#include <esf/sf_uncertainty.hpp>

namespace {

constexpr double kPi = 3.14159265358979323846;


bool close(
    double a,
    double b,
    double rtol = 1e-9,
    double atol = 1e-9
)
{
    return std::abs(a - b) <=
        atol + rtol * std::abs(b);
}


void check_close(double actual, double expected)
{
    assert(close(actual, expected));
}


double mean_of(const std::vector<double>& values)
{
    double sum = 0.0;

    for (const double value : values) {
        sum += value;
    }

    return sum / static_cast<double>(values.size());
}


double sample_var(const std::vector<double>& values)
{
    const double mean = mean_of(values);

    double sum_squared_deviation = 0.0;

    for (const double value : values) {
        const double deviation = value - mean;
        sum_squared_deviation += deviation * deviation;
    }

    return sum_squared_deviation /
        static_cast<double>(values.size() - 1);
}


void test_off_by_default()
{
    agnsf::LightCurve data(
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 1.0, 1.0, 2.0},
        {0.2, 0.2, 0.2, 0.2}
    );

    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(data, bins);

    for (const auto& bin : result.bins()) {
        assert(!bin.measurement.estimated());
        assert(!bin.sampling.estimated());
    }
}


void test_second_order()
{
    agnsf::LightCurve data(
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 1.0, 1.0, 2.0},
        {0.2, 0.2, 0.2, 0.2}
    );

    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::esf::UncertaintyConfig config;
    config.measurement = agnsf::esf::UncertaintyMethod::Analytic;

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            config
        );

    /*
     * Bin 0 (lag 1): pairs (0,1), (1,2), (2,3).
     *
     *   delta = 1, 0, 1; noise = 0.08 each
     *   x_k   = 0.92, -0.08, 0.92
     *
     * mean = 0.586666..., sample_var = 1/3, se = 1/3.
     * sf^2 = 0.586666...; interval on sf via sqrt.
     */
    {
        const std::vector<double> x = {0.92, -0.08, 0.92};

        const double mean = mean_of(x);
        const double se =
            std::sqrt(sample_var(x) / 3.0);

        const double lower =
            std::sqrt(std::max(mean - se, 0.0));
        const double upper =
            std::sqrt(mean + se);

        const auto& bin = result.bin(0);

        assert(bin.measurement.estimated());
        check_close(bin.measurement.lower, lower);
        check_close(bin.measurement.upper, upper);

        // Non-degenerate se gives an asymmetric interval around sf:
        // sqrt is concave, so the two sides differ.
        const double sf = std::sqrt(mean);
        assert(
            std::abs(
                (bin.measurement.upper - sf) -
                (sf - bin.measurement.lower)
            ) > 1e-12
        );
    }

    /*
     * Bin 1 (lag 2): pairs (0,2), (1,3); delta = 1, 1.
     * The true sample variance is zero; a tiny interval width is
     * expected from floating-point cancellation, so a loose
     * tolerance is used here.
     */
    {
        const auto& bin = result.bin(1);

        assert(bin.measurement.estimated());
        assert(close(bin.measurement.lower, bin.measurement.upper, 1e-6));
        assert(close(bin.measurement.lower, std::sqrt(0.92), 1e-6));
    }
}


void test_second_order_no_noise()
{
    agnsf::LightCurve data(
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 1.0, 1.0, 2.0},
        {0.2, 0.2, 0.2, 0.2}
    );

    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::esf::UncertaintyConfig config;
    config.measurement = agnsf::esf::UncertaintyMethod::Analytic;

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrderNoNoise,
            config
        );

    /*
     * Bin 0: x_k = delta^2 = 1, 0, 1.
     *   mean = 2/3, sample_var = 1/3, se = 1/3.
     */
    const std::vector<double> x = {1.0, 0.0, 1.0};

    const double mean = mean_of(x);
    const double se =
        std::sqrt(sample_var(x) / 3.0);

    const auto& bin = result.bin(0);

    check_close(bin.measurement.lower,
        std::sqrt(std::max(mean - se, 0.0)));
    check_close(bin.measurement.upper,
        std::sqrt(mean + se));
}


void test_mean_absolute_deviation()
{
    agnsf::LightCurve data(
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 1.0, 1.0, 2.0},
        {0.2, 0.2, 0.2, 0.2}
    );

    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::esf::UncertaintyConfig config;
    config.measurement = agnsf::esf::UncertaintyMethod::Analytic;

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::MeanAbsoluteDeviation,
            config
        );

    /*
     * Bin 0: |delta| = 1, 0, 1; noise = 0.08 each.
     *
     *   mean_abs = 2/3, sample_var = 1/3, se_abs = 1/3.
     *   d(SF^2) = pi * mean_abs * se_abs = 2*pi/9.
     *   SF^2    = pi/2 * (2/3)^2 - 0.08 = 2*pi/9 - 0.08.
     */
    const double mean_abs = 2.0 / 3.0;
    const double se_abs =
        std::sqrt(sample_var({1.0, 0.0, 1.0}) / 3.0);

    const double delta_sf2 = kPi * mean_abs * se_abs;
    const double sf2 = (kPi / 2.0) * mean_abs * mean_abs - 0.08;

    const auto& bin = result.bin(0);

    check_close(
        bin.measurement.lower,
        std::sqrt(std::max(sf2 - delta_sf2, 0.0))
    );
    check_close(
        bin.measurement.upper,
        std::sqrt(sf2 + delta_sf2)
    );
}


void test_mean_absolute_deviation_no_noise()
{
    agnsf::LightCurve data(
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 1.0, 1.0, 2.0},
        {0.2, 0.2, 0.2, 0.2}
    );

    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::esf::UncertaintyConfig config;
    config.measurement = agnsf::esf::UncertaintyMethod::Analytic;

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::MeanAbsoluteDeviationNoNoise,
            config
        );

    const double mean_abs = 2.0 / 3.0;
    const double se_abs =
        std::sqrt(sample_var({1.0, 0.0, 1.0}) / 3.0);

    const double delta_sf2 = kPi * mean_abs * se_abs;
    const double sf2 = (kPi / 2.0) * mean_abs * mean_abs;

    const auto& bin = result.bin(0);

    check_close(
        bin.measurement.lower,
        std::sqrt(std::max(sf2 - delta_sf2, 0.0))
    );
    check_close(
        bin.measurement.upper,
        std::sqrt(sf2 + delta_sf2)
    );
}


void test_single_pair_bin_unestimated()
{
    /*
     * Only one pair overall -> N = 1 in the bin: the within-bin
     * standard error is undefined, so the interval stays NaN.
     */
    agnsf::LightCurve data(
        {0.0, 3.0},
        {0.0, 1.0},
        {0.1, 0.1}
    );

    const agnsf::esf::LagBins bins({0.0, 4.0});

    agnsf::esf::UncertaintyConfig config;
    config.measurement = agnsf::esf::UncertaintyMethod::Analytic;

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(data, bins, agnsf::esf::SFMethod::SecondOrder, config);

    assert(!result.bin(0).measurement.estimated());
}


void test_noise_dominated_unestimated()
{
    /*
     * noise >> signal: SF^2 < 0, so sf is NaN and no measurement
     * interval can be attached to it.
     */
    agnsf::LightCurve data(
        {0.0, 1.0, 2.0},
        {0.0, 0.0, 0.0},
        {1.0, 1.0, 1.0}
    );

    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::esf::UncertaintyConfig config;
    config.measurement = agnsf::esf::UncertaintyMethod::Analytic;

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            config
        );

    assert(!std::isfinite(result.bin(0).sf));
    assert(!result.bin(0).measurement.estimated());
}


void test_invalid_config()
{
    agnsf::LightCurve data(
        {0.0, 1.0},
        {0.0, 1.0},
        {0.1, 0.1}
    );

    const agnsf::esf::LagBins bins({0.0, 2.0});

    agnsf::esf::SFCalculator calculator;

    // Measurement supports only Off / Analytic.
    {
        agnsf::esf::UncertaintyConfig config;
        config.measurement = agnsf::esf::UncertaintyMethod::Jackknife;

        bool thrown = false;

        try {
            calculator.calculate(data, bins, agnsf::esf::SFMethod::SecondOrder, config);
        }
        catch (const std::invalid_argument&) {
            thrown = true;
        }

        assert(thrown);
    }

    // Sampling is not defined for a single light curve.
    {
        agnsf::esf::UncertaintyConfig config;
        config.sampling = agnsf::esf::UncertaintyMethod::Analytic;

        bool thrown = false;

        try {
            calculator.calculate(data, bins, agnsf::esf::SFMethod::SecondOrder, config);
        }
        catch (const std::invalid_argument&) {
            thrown = true;
        }

        assert(thrown);
    }
}

} // namespace


int main()
{
    test_off_by_default();
    test_second_order();
    test_second_order_no_noise();
    test_mean_absolute_deviation();
    test_mean_absolute_deviation_no_noise();
    test_single_pair_bin_unestimated();
    test_noise_dominated_unestimated();
    test_invalid_config();

    return 0;
}
