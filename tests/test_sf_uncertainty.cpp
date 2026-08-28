#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <core/light_curve.hpp>
#include <esf/sf_calculator.hpp>
#include <core/uncertainty.hpp>

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


agnsf::LightCurve make_curve()
{
    // time {0,1,2,3}, value {0,1,1,2}, error 0.2
    return agnsf::LightCurve(
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 1.0, 1.0, 2.0},
        {0.2, 0.2, 0.2, 0.2}
    );
}


void test_off_by_default()
{
    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(make_curve(), bins);

    for (const auto& bin : result.bins()) {
        assert(!bin.measurement.estimated());
        assert(!bin.within.estimated());
        assert(!bin.sampling.estimated());
    }
}


void test_measurement_second_order()
{
    /*
     * Measurement uncertainty = propagation of the per-observation
     * errors sigma_i (independent-pairs approximation):
     *
     *   Var_meas(SF^2) = 1/N^2 * sum_p [ 4 D_p^2 sigma_{Delta,p}^2
     *                                     + 2 sigma_{Delta,p}^4 ]
     *
     * Bin 0 (lag 1): pairs (0,1),(1,2),(2,3); delta = 1,0,1;
     * noise = sigma_i^2 + sigma_j^2 = 0.08 each.
     *
     *   sum(delta^2 * noise) = 0.16
     *   sum(noise^2)         = 3 * 0.08^2 = 0.0192
     *   Var_meas(SF^2)       = (4*0.16 + 2*0.0192) / 9
     *   sigma_SF             = sigma_SF2 / (2*SF)
     */
    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::UncertaintyConfig config;
    config.measurement = agnsf::UncertaintyMethod::Analytic;

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(
            make_curve(),
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            config
        );

    const double sum_d2_noise = 1.0 * 0.08 + 0.0 * 0.08 + 1.0 * 0.08;
    const double sum_noise2 = 3.0 * 0.08 * 0.08;
    const double var_sf2 =
        (4.0 * sum_d2_noise + 2.0 * sum_noise2) / 9.0;
    const double sigma_sf2 = std::sqrt(var_sf2);

    const double sf2 = (1.0 + 0.0 + 1.0 - 3.0 * 0.08) / 3.0;
    const double sf = std::sqrt(sf2);

    const double sigma_sf = sigma_sf2 / (2.0 * sf);

    const auto& bin = result.bin(0);

    assert(bin.measurement.estimated());
    check_close(bin.measurement.lower, sf - sigma_sf);
    check_close(bin.measurement.upper, sf + sigma_sf);
}


void test_measurement_mean_absolute_deviation()
{
    /*
     * Linear propagation through the mean absolute difference:
     *
     *   Var_meas(SF^2) ~= pi^2 * <|delta|>^2 * sum(sigma_{Delta}^2) / N^2
     *
     * Bin 0: |delta| = 1,0,1; noise = 0.08 each.
     *
     *   mean_abs = 2/3, sum(noise) = 0.24
     */
    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::UncertaintyConfig config;
    config.measurement = agnsf::UncertaintyMethod::Analytic;

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(
            make_curve(),
            bins,
            agnsf::esf::SFMethod::MeanAbsoluteDeviation,
            config
        );

    const double mean_abs = 2.0 / 3.0;
    const double var_sf2 =
        kPi * kPi * mean_abs * mean_abs * 0.24 / 9.0;
    const double sigma_sf2 = std::sqrt(var_sf2);

    const double sf2 = (kPi / 2.0) * mean_abs * mean_abs - 0.08;
    const double sf = std::sqrt(sf2);

    const double sigma_sf = sigma_sf2 / (2.0 * sf);

    const auto& bin = result.bin(0);

    assert(bin.measurement.estimated());
    check_close(bin.measurement.lower, sf - sigma_sf);
    check_close(bin.measurement.upper, sf + sigma_sf);
}


void test_measurement_monte_carlo()
{
    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::UncertaintyConfig config;
    config.measurement = agnsf::UncertaintyMethod::MonteCarlo;
    config.n_bootstrap = 500;
    config.bootstrap_seed = 42;

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult first =
        calculator.calculate(
            make_curve(),
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            config
        );

    const agnsf::esf::SFResult second =
        calculator.calculate(
            make_curve(),
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            config
        );

    // Reproducible with the same seed.
    check_close(first.bin(0).measurement.lower,
                second.bin(0).measurement.lower);
    check_close(first.bin(0).measurement.upper,
                second.bin(0).measurement.upper);

    const auto& bin = first.bin(0);

    assert(bin.measurement.estimated());
    assert(std::isfinite(bin.measurement.lower));
    assert(std::isfinite(bin.measurement.upper));
    assert(bin.measurement.lower <= bin.measurement.upper);
    assert(bin.measurement.lower < bin.measurement.upper);  // noise -> width > 0

    // The point estimate should lie within the 16/84 interval for
    // this symmetric-ish data.
    assert(bin.measurement.lower <= bin.sf + 1e-9);
    assert(bin.measurement.upper >= bin.sf - 1e-9);
}


void test_within_second_order()
{
    /*
     * Naive within-bin statistical uncertainty:
     *
     *   x_k = delta_k^2 - noise_k;  se = sample_std(x) / sqrt(N)
     *
     * Bin 0: x = 0.92, -0.08, 0.92; mean = 0.5866667,
     * sample_var = 1/3, se = 1/3.
     */
    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::UncertaintyConfig config;
    config.within = agnsf::UncertaintyMethod::Analytic;

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(
            make_curve(),
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            config
        );

    const std::vector<double> x = {0.92, -0.08, 0.92};

    const double mean = (0.92 - 0.08 + 0.92) / 3.0;
    const double sample_var =
        (
            (0.92 - mean) * (0.92 - mean) +
            (-0.08 - mean) * (-0.08 - mean) +
            (0.92 - mean) * (0.92 - mean)
        ) / 2.0;
    const double se = std::sqrt(sample_var / 3.0);

    const auto& bin = result.bin(0);

    assert(bin.within.estimated());
    check_close(
        bin.within.lower,
        std::sqrt(std::max(mean - se, 0.0))
    );
    check_close(
        bin.within.upper,
        std::sqrt(mean + se)
    );
}


void test_within_mean_absolute_deviation()
{
    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::UncertaintyConfig config;
    config.within = agnsf::UncertaintyMethod::Analytic;

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(
            make_curve(),
            bins,
            agnsf::esf::SFMethod::MeanAbsoluteDeviation,
            config
        );

    const double mean_abs = 2.0 / 3.0;
    const double sample_var =
        (
            (1.0 - mean_abs) * (1.0 - mean_abs) +
            (0.0 - mean_abs) * (0.0 - mean_abs) +
            (1.0 - mean_abs) * (1.0 - mean_abs)
        ) / 2.0;
    const double se_abs = std::sqrt(sample_var / 3.0);
    const double delta_sf2 = kPi * mean_abs * se_abs;

    const double sf2 = (kPi / 2.0) * mean_abs * mean_abs - 0.08;

    const auto& bin = result.bin(0);

    assert(bin.within.estimated());
    check_close(
        bin.within.lower,
        std::sqrt(std::max(sf2 - delta_sf2, 0.0))
    );
    check_close(
        bin.within.upper,
        std::sqrt(sf2 + delta_sf2)
    );
}


void test_noise_dominated_unestimated()
{
    // noise >> signal -> SF^2 < 0 -> sf = NaN -> no measurement.
    agnsf::LightCurve data(
        {0.0, 1.0, 2.0},
        {0.0, 0.0, 0.0},
        {1.0, 1.0, 1.0}
    );

    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::UncertaintyConfig config;
    config.measurement = agnsf::UncertaintyMethod::Analytic;

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
    const agnsf::esf::LagBins bins({0.0, 2.0});

    agnsf::esf::SFCalculator calculator;

    const auto expect_throw =
        [&](const agnsf::UncertaintyConfig& config)
        {
            bool thrown = false;

            try {
                calculator.calculate(
                    make_curve(),
                    bins,
                    agnsf::esf::SFMethod::SecondOrder,
                    config
                );
            }
            catch (const std::invalid_argument&) {
                thrown = true;
            }

            assert(thrown);
        };

    // measurement does not support Jackknife / Bootstrap.
    {
        agnsf::UncertaintyConfig config;
        config.measurement = agnsf::UncertaintyMethod::Jackknife;
        expect_throw(config);
    }

    // within does not support Jackknife / Bootstrap / MonteCarlo.
    {
        agnsf::UncertaintyConfig config;
        config.within = agnsf::UncertaintyMethod::Jackknife;
        expect_throw(config);
    }

    // sampling is not defined for a single light curve.
    {
        agnsf::UncertaintyConfig config;
        config.sampling = agnsf::UncertaintyMethod::Analytic;
        expect_throw(config);
    }
}

} // namespace


int main()
{
    test_off_by_default();
    test_measurement_second_order();
    test_measurement_mean_absolute_deviation();
    test_measurement_monte_carlo();
    test_within_second_order();
    test_within_mean_absolute_deviation();
    test_noise_dominated_unestimated();
    test_invalid_config();

    return 0;
}
