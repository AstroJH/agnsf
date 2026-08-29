#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

#include <core/light_curve.hpp>
#include <timedelay/cross_correlation.hpp>
#include <timedelay/fr_rss.hpp>

namespace {

const double kPi = 3.14159265358979323846;


bool close(
    double a,
    double b,
    double atol = 1e-9
)
{
    return std::abs(a - b) <= atol;
}


/**
 * A delayed pair of light curves on a common regular grid:
 *
 *   f(t) = 3 sin(2 pi t / 50) + 1.5 sin(2 pi t / 23)
 *   g(t) = f(t - lag_true) + noise
 */
void make_delayed_pair(
    agnsf::LightCurve& continuum,
    agnsf::LightCurve& response,
    double lag_true = 7.0,
    double noise = 0.05,
    int n = 120
)
{
    std::vector<double> t;
    std::vector<double> f;
    std::vector<double> g;
    std::vector<double> e;

    const auto signal = [](double x) {
        return 3.0 * std::sin(2.0 * kPi * x / 50.0) +
               1.5 * std::sin(2.0 * kPi * x / 23.0);
    };

    // Simple deterministic pseudo-noise.
    std::uint32_t state = 12345;

    const auto rnd = [&state]() {
        state = state * 1664525u + 1013904223u;
        return static_cast<double>((state >> 8) & 0xffffffu) /
               static_cast<double>(0x1000000u) - 0.5;
    };

    for (int i = 0; i < n; ++i) {
        const double ti = static_cast<double>(i);
        const double fi = signal(ti);
        const double gi = signal(ti - lag_true) + noise * rnd();

        t.push_back(ti);
        f.push_back(fi);
        g.push_back(gi);
        e.push_back(noise);
    }

    continuum = agnsf::LightCurve(t, f, e);
    response = agnsf::LightCurve(t, g, e);
}


void test_dcf_peak()
{
    agnsf::LightCurve c({}, {}, {});
    agnsf::LightCurve r({}, {}, {});

    make_delayed_pair(c, r, /*lag_true=*/7.0);

    agnsf::timedelay::LagGrid grid;
    grid.min = -30.0;
    grid.max = 30.0;
    grid.step = 1.0;

    agnsf::timedelay::CrossCorrelationConfig config;
    config.method = agnsf::timedelay::CrossCorrelationMethod::Dcf;

    const auto result =
        agnsf::timedelay::cross_correlate(c, r, grid, config);

    assert(result.tau.size() == 61);
    assert(std::isfinite(result.lag_peak));
    assert(std::abs(result.lag_peak - 7.0) <= 1.0);
    assert(result.peak_value > 0.9);
}


void test_iccf_peak_and_centroid()
{
    agnsf::LightCurve c({}, {}, {});
    agnsf::LightCurve r({}, {}, {});

    make_delayed_pair(c, r, /*lag_true=*/7.0);

    agnsf::timedelay::LagGrid grid;
    grid.min = -30.0;
    grid.max = 30.0;
    grid.step = 1.0;

    agnsf::timedelay::CrossCorrelationConfig config;
    config.method = agnsf::timedelay::CrossCorrelationMethod::Iccf;

    const auto result =
        agnsf::timedelay::cross_correlate(c, r, grid, config);

    assert(std::isfinite(result.lag_peak));
    assert(std::abs(result.lag_peak - 7.0) <= 1.0);
    assert(result.peak_value > 0.9);

    // The centroid should also be near the true lag.
    assert(std::isfinite(result.lag_centroid));
    assert(std::abs(result.lag_centroid - 7.0) <= 2.0);

    // lag_value selects the requested estimate.
    assert(
        close(
            agnsf::timedelay::lag_value(
                result,
                agnsf::timedelay::LagEstimate::Peak
            ),
            result.lag_peak
        )
    );

    assert(
        close(
            agnsf::timedelay::lag_value(
                result,
                agnsf::timedelay::LagEstimate::Centroid
            ),
            result.lag_centroid
        )
    );
}


void test_fr_rss()
{
    agnsf::LightCurve c({}, {}, {});
    agnsf::LightCurve r({}, {}, {});

    make_delayed_pair(c, r, /*lag_true=*/7.0, /*noise=*/0.05);

    agnsf::timedelay::LagGrid grid;
    grid.min = -30.0;
    grid.max = 30.0;
    grid.step = 1.0;

    agnsf::timedelay::CrossCorrelationConfig ccf_config;
    ccf_config.method = agnsf::timedelay::CrossCorrelationMethod::Iccf;

    agnsf::timedelay::FRRSSConfig config;
    config.n_realizations = 300;
    config.seed = 42;

    const agnsf::Uncertainty first =
        agnsf::timedelay::lag_uncertainty(
            c,
            r,
            grid,
            agnsf::timedelay::LagEstimate::Peak,
            ccf_config,
            config
        );

    const agnsf::Uncertainty second =
        agnsf::timedelay::lag_uncertainty(
            c,
            r,
            grid,
            agnsf::timedelay::LagEstimate::Peak,
            ccf_config,
            config
        );

    // Reproducible with the same seed.
    assert(close(first.lower, second.lower));
    assert(close(first.upper, second.upper));

    assert(first.estimated());
    assert(std::isfinite(first.lower));
    assert(std::isfinite(first.upper));
    assert(first.lower <= first.upper);

    // For a strong signal the interval brackets the true lag.
    assert(first.lower <= 8.0);
    assert(first.upper >= 6.0);
}


void test_invalid_inputs()
{
    agnsf::LightCurve c({}, {}, {});
    agnsf::LightCurve r({}, {}, {});

    make_delayed_pair(c, r, 7.0);

    // Invalid grid.
    agnsf::timedelay::LagGrid bad_grid;
    bad_grid.min = 10.0;
    bad_grid.max = 5.0;

    bool thrown = false;

    try {
        agnsf::timedelay::cross_correlate(c, r, bad_grid);
    }
    catch (const std::invalid_argument&) {
        thrown = true;
    }

    assert(thrown);

    // Empty curve.
    agnsf::LightCurve empty({}, {}, {});

    thrown = false;

    try {
        agnsf::timedelay::cross_correlate(
            empty,
            r,
            agnsf::timedelay::LagGrid{}
        );
    }
    catch (const std::invalid_argument&) {
        thrown = true;
    }

    assert(thrown);
}

} // namespace


int main()
{
    test_dcf_peak();
    test_iccf_peak_and_centroid();
    test_fr_rss();
    test_invalid_inputs();

    return 0;
}
