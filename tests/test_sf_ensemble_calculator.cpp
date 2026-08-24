#include <cassert>
#include <cmath>
#include <vector>

#include <esf/lag_bins.hpp>
#include <esf/light_curve.hpp>
#include <esf/sf_ensemble_calculator.hpp>

namespace {

void test_mean_of_individual_sfs()
{
    esf::LightCurve lc1(
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 0.0}
    );

    esf::LightCurve lc2(
        {0.0, 1.0},
        {0.0, 2.0},
        {0.0, 0.0}
    );

    std::vector<esf::LightCurve> data{
        lc1,
        lc2
    };

    esf::LagBins bins({
        0.0,
        2.0
    });

    esf::SFEnsembleCalculator calculator;

    const auto result =
        calculator.calculate(data, bins);

    assert(result.size() == 1);

    const auto& bin = result.bin(0);

    /*
     * Individual SFs:
     *
     * LC1: SF = 1
     * LC2: SF = 2
     *
     * Mean:
     *
     * ESF = (1 + 2) / 2 = 1.5
     */
    assert(bin.count == 2);

    assert(
        std::abs(bin.sf - 1.5)
        < 1e-12
    );

    assert(
        std::abs(bin.sf_squared - 2.25)
        < 1e-12
    );
}

void test_pooled_and_ensemble_are_different()
{
    esf::LightCurve lc1(
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 0.0}
    );

    esf::LightCurve lc2(
        {0.0, 1.0},
        {0.0, 2.0},
        {0.0, 0.0}
    );

    std::vector<esf::LightCurve> data{
        lc1,
        lc2
    };

    esf::LagBins bins({
        0.0,
        2.0
    });

    esf::SFEnsembleCalculator calculator;

    const auto result =
        calculator.calculate(data, bins);

    const auto& bin = result.bin(0);

    // Arithmetic mean of individual SFs:
    //
    // (1 + 2) / 2 = 1.5
    assert(
        std::abs(bin.sf - 1.5) < 1e-12
    );

    // This is intentionally different from:
    //
    // sqrt((1^2 + 2^2) / 2) = sqrt(2.5)
    assert(
        std::abs(
            bin.sf - std::sqrt(2.5)
        ) > 1e-6
    );
}
} // namespace


int main()
{
    test_mean_of_individual_sfs();
    test_pooled_and_ensemble_are_different();

    return 0;
}