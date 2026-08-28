#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <core/light_curve.hpp>
#include <esf/pooled_esf_calculator.hpp>
#include <esf/sf_calculator.hpp>
#include <esf/sf_ensemble_calculator.hpp>
#include <esf/sf_uncertainty.hpp>

namespace {

bool close(
    double a,
    double b,
    double rtol = 1e-8,
    double atol = 1e-8
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


double sample_std(const std::vector<double>& values)
{
    const double mean = mean_of(values);

    double sum_squared_deviation = 0.0;

    for (const double value : values) {
        const double deviation = value - mean;
        sum_squared_deviation += deviation * deviation;
    }

    return std::sqrt(
        sum_squared_deviation /
        static_cast<double>(values.size() - 1)
    );
}


/**
 * Independent jackknife interval over the given values.
 * `sqrt_mean` mirrors the SqrtMeanSquared final transform.
 */
void jackknife_expected(
    const std::vector<double>& values,
    bool sqrt_mean,
    double& lower,
    double& upper
)
{
    const double n =
        static_cast<double>(values.size());

    double total = 0.0;

    for (const double value : values) {
        total += value;
    }

    std::vector<double> leave_one_out;

    for (const double value : values) {
        const double mean_without =
            (total - value) / (n - 1.0);

        const double statistic =
            sqrt_mean
                ? std::sqrt(std::max(mean_without, 0.0))
                : mean_without;

        leave_one_out.push_back(statistic);
    }

    const double leave_one_out_mean =
        mean_of(leave_one_out);

    double sum_squared_deviation = 0.0;

    for (const double value : leave_one_out) {
        const double deviation =
            value - leave_one_out_mean;

        sum_squared_deviation +=
            deviation * deviation;
    }

    const double variance =
        ((n - 1.0) / n) * sum_squared_deviation;

    const double sigma = std::sqrt(variance);

    const double full_mean = total / n;

    const double point =
        sqrt_mean
            ? std::sqrt(std::max(full_mean, 0.0))
            : full_mean;

    lower = point - sigma;
    upper = point + sigma;
}


std::vector<agnsf::LightCurve> make_curves()
{
    // Three curves; per-curve bin0 (lag 1) SF^2 = 0.98 / 3.98 / 8.98.
    return {
        agnsf::LightCurve(
            {0.0, 1.0, 2.0},
            {0.0, 1.0, 2.0},
            {0.1, 0.1, 0.1}
        ),
        agnsf::LightCurve(
            {0.0, 1.0, 2.0},
            {0.0, 2.0, 4.0},
            {0.1, 0.1, 0.1}
        ),
        agnsf::LightCurve(
            {0.0, 1.0, 2.0},
            {0.0, 3.0, 6.0},
            {0.1, 0.1, 0.1}
        )
    };
}


void test_analytic_sampling_sqrt_mean_squared()
{
    const std::vector<agnsf::LightCurve> data =
        make_curves();

    const agnsf::esf::LagBins bins({0.0, 1.5});

    agnsf::esf::UncertaintyConfig config;
    config.sampling = agnsf::esf::UncertaintyMethod::Analytic;

    agnsf::esf::SFEnsembleCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            agnsf::esf::SFEnsembleCalculator::Method::SqrtMeanSquared,
            config
        );

    /*
     * Per-curve SF^2 = {0.98, 3.98, 8.98}; mean = 4.646667.
     * se = sample_std / sqrt(3); interval on sf^2 mapped through sqrt.
     */
    const std::vector<double> values = {0.98, 3.98, 8.98};

    const double mean = mean_of(values);
    const double se =
        sample_std(values) / std::sqrt(3.0);

    const auto& bin = result.bin(0);

    assert(bin.sampling.estimated());
    check_close(
        bin.sampling.lower,
        std::sqrt(std::max(mean - se, 0.0))
    );
    check_close(
        bin.sampling.upper,
        std::sqrt(mean + se)
    );
}


void test_analytic_sampling_mean_sf()
{
    const std::vector<agnsf::LightCurve> data =
        make_curves();

    const agnsf::esf::LagBins bins({0.0, 1.5});

    agnsf::esf::UncertaintyConfig config;
    config.sampling = agnsf::esf::UncertaintyMethod::Analytic;

    agnsf::esf::SFCalculator sf_calculator;

    const std::vector<double> sf_values = {
        sf_calculator.calculate(data[0], bins).bin(0).sf,
        sf_calculator.calculate(data[1], bins).bin(0).sf,
        sf_calculator.calculate(data[2], bins).bin(0).sf
    };

    agnsf::esf::SFEnsembleCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            agnsf::esf::SFEnsembleCalculator::Method::MeanSf,
            config
        );

    const double mean = mean_of(sf_values);
    const double se =
        sample_std(sf_values) / std::sqrt(3.0);

    const auto& bin = result.bin(0);

    check_close(bin.sampling.lower, mean - se);
    check_close(bin.sampling.upper, mean + se);
}


void test_jackknife_sampling()
{
    const std::vector<agnsf::LightCurve> data =
        make_curves();

    const agnsf::esf::LagBins bins({0.0, 1.5});

    agnsf::esf::UncertaintyConfig config;
    config.sampling = agnsf::esf::UncertaintyMethod::Jackknife;

    agnsf::esf::SFEnsembleCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            agnsf::esf::SFEnsembleCalculator::Method::SqrtMeanSquared,
            config
        );

    double expected_lower = 0.0;
    double expected_upper = 0.0;

    jackknife_expected(
        {0.98, 3.98, 8.98},
        true,
        expected_lower,
        expected_upper
    );

    const auto& bin = result.bin(0);

    check_close(bin.sampling.lower, expected_lower);
    check_close(bin.sampling.upper, expected_upper);
}


void test_bootstrap_sampling_reproducible()
{
    const std::vector<agnsf::LightCurve> data =
        make_curves();

    const agnsf::esf::LagBins bins({0.0, 1.5});

    agnsf::esf::UncertaintyConfig config;
    config.sampling = agnsf::esf::UncertaintyMethod::Bootstrap;
    config.n_bootstrap = 200;
    config.bootstrap_seed = 12345;

    agnsf::esf::SFEnsembleCalculator calculator;

    const agnsf::esf::SFResult first =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            agnsf::esf::SFEnsembleCalculator::Method::MeanSf,
            config
        );

    const agnsf::esf::SFResult second =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            agnsf::esf::SFEnsembleCalculator::Method::MeanSf,
            config
        );

    const auto& a = first.bin(0).sampling;
    const auto& b = second.bin(0).sampling;

    assert(a.estimated());
    check_close(a.lower, b.lower);
    check_close(a.upper, b.upper);

    // The interval must bracket the point estimate.
    const double sf = first.bin(0).sf;

    assert(a.lower <= sf + 1e-9);
    assert(a.upper >= sf - 1e-9);
}


void test_measurement_propagation_mean_sf()
{
    // Curves with in-bin scatter so the per-curve measurement is
    // non-degenerate.
    std::vector<agnsf::LightCurve> data = {
        agnsf::LightCurve(
            {0.0, 1.0, 2.0, 3.0},
            {0.0, 1.0, 1.0, 2.0},
            {0.2, 0.2, 0.2, 0.2}
        ),
        agnsf::LightCurve(
            {0.0, 1.0, 2.0, 3.0},
            {0.0, 2.0, 1.0, 3.0},
            {0.2, 0.2, 0.2, 0.2}
        ),
        agnsf::LightCurve(
            {0.0, 1.0, 2.0, 3.0},
            {0.0, 3.0, 2.0, 4.0},
            {0.2, 0.2, 0.2, 0.2}
        )
    };

    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::esf::UncertaintyConfig per_curve_config;
    per_curve_config.measurement =
        agnsf::esf::UncertaintyMethod::Analytic;

    agnsf::esf::SFCalculator sf_calculator;

    // Per-curve measurement half-widths (MeanSf combines them in
    // quadrature and divides by the number of curves).
    std::vector<double> sf_values;
    std::vector<double> half_widths;

    for (const auto& curve : data) {

        const agnsf::esf::SFResult r =
            sf_calculator.calculate(
                curve,
                bins,
                agnsf::esf::SFMethod::SecondOrder,
                per_curve_config
            );

        const auto& measurement = r.bin(0).measurement;

        assert(measurement.estimated());

        sf_values.push_back(r.bin(0).sf);

        half_widths.push_back(
            (measurement.upper - measurement.lower) / 2.0
        );
    }

    const double n =
        static_cast<double>(sf_values.size());

    const double esf = mean_of(sf_values);

    double sum_squared_sigma = 0.0;

    for (const double half_width : half_widths) {
        sum_squared_sigma += half_width * half_width;
    }

    const double sigma =
        std::sqrt(sum_squared_sigma) / n;

    agnsf::esf::UncertaintyConfig config;
    config.measurement = agnsf::esf::UncertaintyMethod::Analytic;

    agnsf::esf::SFEnsembleCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            agnsf::esf::SFEnsembleCalculator::Method::MeanSf,
            config
        );

    const auto& measurement = result.bin(0).measurement;

    assert(measurement.estimated());
    check_close(measurement.lower, esf - sigma);
    check_close(measurement.upper, esf + sigma);
}


void test_pooled_jackknife()
{
    /*
     * Three single-pair curves; pooled SF^2 = (0.98+3.98+8.98)/3.
     * Leave-one-out pooled SFs are computed independently below.
     */
    std::vector<agnsf::LightCurve> data = {
        agnsf::LightCurve({0.0, 1.0}, {0.0, 1.0}, {0.1, 0.1}),
        agnsf::LightCurve({0.0, 1.0}, {0.0, 2.0}, {0.1, 0.1}),
        agnsf::LightCurve({0.0, 1.0}, {0.0, 3.0}, {0.1, 0.1})
    };

    const agnsf::esf::LagBins bins({0.0, 1.5});

    agnsf::esf::UncertaintyConfig config;
    config.sampling = agnsf::esf::UncertaintyMethod::Jackknife;

    agnsf::esf::PooledESFCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(data, bins, agnsf::esf::SFMethod::SecondOrder, config);

    const auto& bin = result.bin(0);

    // Independent leave-one-out pooled sf values.
    const std::vector<double> loo_sf = {
        std::sqrt((3.98 + 8.98) / 2.0),
        std::sqrt((0.98 + 8.98) / 2.0),
        std::sqrt((0.98 + 3.98) / 2.0)
    };

    double expected_lower = 0.0;
    double expected_upper = 0.0;

    jackknife_expected(loo_sf, false, expected_lower, expected_upper);

    check_close(bin.sampling.lower, expected_lower);
    check_close(bin.sampling.upper, expected_upper);

    // The point estimate must be the pooled sf of all three curves.
    check_close(bin.sf, std::sqrt(13.94 / 3.0));
}


void test_pooled_bootstrap_reproducible()
{
    std::vector<agnsf::LightCurve> data = {
        agnsf::LightCurve({0.0, 1.0}, {0.0, 1.0}, {0.1, 0.1}),
        agnsf::LightCurve({0.0, 1.0}, {0.0, 2.0}, {0.1, 0.1}),
        agnsf::LightCurve({0.0, 1.0}, {0.0, 3.0}, {0.1, 0.1})
    };

    const agnsf::esf::LagBins bins({0.0, 1.5});

    agnsf::esf::UncertaintyConfig config;
    config.sampling = agnsf::esf::UncertaintyMethod::Bootstrap;
    config.n_bootstrap = 300;
    config.bootstrap_seed = 7;

    agnsf::esf::PooledESFCalculator calculator;

    const agnsf::esf::SFResult first =
        calculator.calculate(data, bins, agnsf::esf::SFMethod::SecondOrder, config);

    const agnsf::esf::SFResult second =
        calculator.calculate(data, bins, agnsf::esf::SFMethod::SecondOrder, config);

    const auto& a = first.bin(0).sampling;
    const auto& b = second.bin(0).sampling;

    assert(a.estimated());
    check_close(a.lower, b.lower);
    check_close(a.upper, b.upper);
    assert(a.lower <= a.upper);
}


void test_pooled_invalid_config()
{
    std::vector<agnsf::LightCurve> data = {
        agnsf::LightCurve({0.0, 1.0}, {0.0, 1.0}, {0.1, 0.1}),
        agnsf::LightCurve({0.0, 1.0}, {0.0, 2.0}, {0.1, 0.1})
    };

    const agnsf::esf::LagBins bins({0.0, 1.5});

    agnsf::esf::PooledESFCalculator calculator;

    // Analytic sampling is not defined for pooled ESF.
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

    // Bootstrap needs at least one replicate.
    {
        agnsf::esf::UncertaintyConfig config;
        config.sampling = agnsf::esf::UncertaintyMethod::Bootstrap;
        config.n_bootstrap = 0;

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


void test_within_propagation_mean_sf()
{
    /*
     * The naive within-bin statistical uncertainty of each curve is
     * combined in quadrature across independent curves (same rule as
     * measurement propagation).
     */
    std::vector<agnsf::LightCurve> data = {
        agnsf::LightCurve(
            {0.0, 1.0, 2.0, 3.0},
            {0.0, 1.0, 1.0, 2.0},
            {0.2, 0.2, 0.2, 0.2}
        ),
        agnsf::LightCurve(
            {0.0, 1.0, 2.0, 3.0},
            {0.0, 2.0, 1.0, 3.0},
            {0.2, 0.2, 0.2, 0.2}
        )
    };

    const agnsf::esf::LagBins bins({0.0, 1.5, 3.0});

    agnsf::esf::UncertaintyConfig per_curve_config;
    per_curve_config.within = agnsf::esf::UncertaintyMethod::Analytic;

    agnsf::esf::SFCalculator sf_calculator;

    std::vector<double> sf_values;
    std::vector<double> half_widths;

    for (const auto& curve : data) {

        const agnsf::esf::SFResult r =
            sf_calculator.calculate(
                curve,
                bins,
                agnsf::esf::SFMethod::SecondOrder,
                per_curve_config
            );

        const auto& within = r.bin(0).within;

        assert(within.estimated());

        sf_values.push_back(r.bin(0).sf);

        half_widths.push_back(
            (within.upper - within.lower) / 2.0
        );
    }

    const double n = static_cast<double>(sf_values.size());
    const double esf = (
        sf_values[0] + sf_values[1]
    ) / n;

    double sum_squared_sigma = 0.0;

    for (const double half_width : half_widths) {
        sum_squared_sigma += half_width * half_width;
    }

    const double sigma =
        std::sqrt(sum_squared_sigma) / n;

    agnsf::esf::UncertaintyConfig config;
    config.within = agnsf::esf::UncertaintyMethod::Analytic;

    agnsf::esf::SFEnsembleCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            agnsf::esf::SFEnsembleCalculator::Method::MeanSf,
            config
        );

    const auto& within = result.bin(0).within;

    assert(within.estimated());
    check_close(within.lower, esf - sigma);
    check_close(within.upper, esf + sigma);
}


void test_pooled_measurement_monte_carlo()
{
    std::vector<agnsf::LightCurve> data = {
        agnsf::LightCurve({0.0, 1.0}, {0.0, 1.0}, {0.1, 0.1}),
        agnsf::LightCurve({0.0, 1.0}, {0.0, 2.0}, {0.1, 0.1})
    };

    const agnsf::esf::LagBins bins({0.0, 1.5});

    agnsf::esf::UncertaintyConfig config;
    config.measurement = agnsf::esf::UncertaintyMethod::MonteCarlo;
    config.n_bootstrap = 400;
    config.bootstrap_seed = 7;

    agnsf::esf::PooledESFCalculator calculator;

    const agnsf::esf::SFResult first =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            config
        );

    const agnsf::esf::SFResult second =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            config
        );

    check_close(first.bin(0).measurement.lower,
                second.bin(0).measurement.lower);
    check_close(first.bin(0).measurement.upper,
                second.bin(0).measurement.upper);

    const auto& bin = first.bin(0);

    assert(bin.measurement.estimated());
    assert(std::isfinite(bin.measurement.lower));
    assert(std::isfinite(bin.measurement.upper));
    assert(bin.measurement.lower <= bin.measurement.upper);
}

} // namespace


int main()
{
    test_analytic_sampling_sqrt_mean_squared();
    test_analytic_sampling_mean_sf();
    test_jackknife_sampling();
    test_bootstrap_sampling_reproducible();
    test_measurement_propagation_mean_sf();
    test_pooled_jackknife();
    test_pooled_bootstrap_reproducible();
    test_pooled_invalid_config();
    test_within_propagation_mean_sf();
    test_pooled_measurement_monte_carlo();

    return 0;
}
