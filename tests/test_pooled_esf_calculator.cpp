#include <cassert>
#include <cmath>
#include <vector>

#include <esf/lag_bins.hpp>
#include <core/light_curve.hpp>
#include <esf/pooled_esf_calculator.hpp>

namespace {

void test_pooled_pairs()
{
    agnsf::LightCurve lc1(
        {0.0, 1.0},
        {0.0, 1.0},
        {0.0, 0.0}
    );

    agnsf::LightCurve lc2(
        {0.0, 1.0},
        {0.0, 2.0},
        {0.0, 0.0}
    );

    std::vector<agnsf::LightCurve> data{
        lc1,
        lc2
    };

    agnsf::esf::LagBins bins({
        0.0,
        2.0
    });

    agnsf::esf::PooledESFCalculator calculator;

    const auto result =
        calculator.calculate(data, bins);

    assert(result.size() == 1);

    const auto& bin = result.bin(0);

    /*
     * LC1:
     * Δ² = 1
     *
     * LC2:
     * Δ² = 4
     *
     * pooled SF² = (1 + 4) / 2 = 2.5
     */
    assert(bin.count == 2);

    assert(
        std::abs(bin.sf_squared - 2.5)
        < 1e-12
    );

    assert(
        std::abs(bin.sf - std::sqrt(2.5))
        < 1e-12
    );
}


void test_pooled_methods()
{
    constexpr double kPi = 3.14159265358979323846;

    agnsf::LightCurve lc1(
        {0.0, 1.0},
        {0.0, 1.0},
        {0.1, 0.1}
    );

    agnsf::LightCurve lc2(
        {0.0, 1.0},
        {0.0, 2.0},
        {0.1, 0.1}
    );

    std::vector<agnsf::LightCurve> data{
        lc1,
        lc2
    };

    agnsf::esf::LagBins bins({
        0.0,
        2.0
    });

    agnsf::esf::PooledESFCalculator calculator;

    /*
     * Pooled statistics:
     *
     *   sum(delta^2) = 1 + 4 = 5
     *   sum(|delta|) = 1 + 2 = 3
     *   sum(noise)   = 0.02 + 0.02 = 0.04
     *   n            = 2
     */
    {
        const auto result =
            calculator.calculate(
                data,
                bins,
                agnsf::esf::SFMethod::SecondOrder
            );

        const auto& bin = result.bin(0);

        assert(bin.count == 2);
        assert(
            std::abs(
                bin.sf_squared - 2.48
            ) < 1e-12
        );
    }

    {
        const auto result =
            calculator.calculate(
                data,
                bins,
                agnsf::esf::SFMethod::SecondOrderNoNoise
            );

        const auto& bin = result.bin(0);

        assert(bin.count == 2);
        assert(
            std::abs(
                bin.sf_squared - 2.5
            ) < 1e-12
        );
    }

    {
        const auto result =
            calculator.calculate(
                data,
                bins,
                agnsf::esf::SFMethod::MeanAbsoluteDeviation
            );

        const auto& bin = result.bin(0);

        // pi/2 * (3/2)^2 - 0.04/2 = pi/2 * 2.25 - 0.02
        assert(bin.count == 2);
        assert(
            std::abs(
                bin.sf_squared -
                (kPi / 2.0 * 2.25 - 0.02)
            ) < 1e-12
        );
    }

    {
        const auto result =
            calculator.calculate(
                data,
                bins,
                agnsf::esf::SFMethod::MeanAbsoluteDeviationNoNoise
            );

        const auto& bin = result.bin(0);

        // pi/2 * (3/2)^2 = pi/2 * 2.25
        assert(bin.count == 2);
        assert(
            std::abs(
                bin.sf_squared -
                (kPi / 2.0 * 2.25)
            ) < 1e-12
        );
    }
}

} // namespace


int main()
{
    test_pooled_pairs();
    test_pooled_methods();

    return 0;
}