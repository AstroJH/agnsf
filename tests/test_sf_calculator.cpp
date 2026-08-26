#include <cassert>
#include <cmath>

#include <esf/lag_bins.hpp>
#include <core/light_curve.hpp>
#include <esf/sf_calculator.hpp>

namespace {

void test_simple_linear_signal()
{
    agnsf::LightCurve data(
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 0.0, 0.0, 0.0}
    );

    agnsf::esf::LagBins bins({
        0.0,
        1.5,
        2.5,
        4.0
    });

    agnsf::esf::SFCalculator calculator;

    const auto result =
        calculator.calculate(data, bins);

    assert(result.size() == 3);

    {
        const auto& bin = result.bin(0);

        assert(bin.count == 3);

        assert(
            std::abs(bin.sf_squared - 1.0)
            < 1e-12
        );

        assert(
            std::abs(bin.sf - 1.0)
            < 1e-12
        );
    }

    {
        const auto& bin = result.bin(1);

        assert(bin.count == 2);

        assert(
            std::abs(bin.sf_squared - 4.0)
            < 1e-12
        );

        assert(
            std::abs(bin.sf - 2.0)
            < 1e-12
        );
    }

    {
        const auto& bin = result.bin(2);

        assert(bin.count == 1);

        assert(
            std::abs(bin.sf_squared - 9.0)
            < 1e-12
        );

        assert(
            std::abs(bin.sf - 3.0)
            < 1e-12
        );
    }
}


void test_noise_correction()
{
    agnsf::LightCurve data(
        {0.0, 1.0, 2.0},
        {0.0, 1.0, 2.0},
        {0.1, 0.1, 0.1}
    );

    agnsf::esf::LagBins bins({
        0.0,
        1.5,
        3.0
    });

    agnsf::esf::SFCalculator calculator;

    const auto result =
        calculator.calculate(data, bins);

    {
        const auto& bin = result.bin(0);

        assert(bin.count == 2);

        assert(
            std::abs(bin.sf_squared - 0.98)
            < 1e-12
        );
    }

    {
        const auto& bin = result.bin(1);

        assert(bin.count == 1);

        assert(
            std::abs(bin.sf_squared - 3.98)
            < 1e-12
        );
    }
}

void test_bin_boundaries()
{
    agnsf::LightCurve data(
        {0.0, 1.0, 2.0},
        {0.0, 1.0, 2.0},
        {0.0, 0.0, 0.0}
    );

    agnsf::esf::LagBins bins({
        0.0,
        1.0,
        2.0,
        3.0
    });

    agnsf::esf::SFCalculator calculator;

    const auto result =
        calculator.calculate(data, bins);

    // tau = 1 belongs to [1, 2)
    assert(result.bin(0).count == 0);
    assert(result.bin(1).count == 2);

    // tau = 2 belongs to [2, 3)
    assert(result.bin(2).count == 1);
}

void test_methods()
{
    constexpr double kPi = 3.14159265358979323846;

    agnsf::LightCurve data(
        {0.0, 1.0, 2.0},
        {0.0, 1.0, 2.0},
        {0.1, 0.1, 0.1}
    );

    agnsf::esf::LagBins bins({
        0.0,
        1.5,
        3.0
    });

    agnsf::esf::SFCalculator calculator;

    const auto result =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::MeanAbsoluteDeviation
        );

    /*
     * Bin 0 (lag 1): pairs (0,1) and (1,2).
     *   sum(delta^2) = 2, sum(|delta|) = 2, sum(noise) = 0.04, n = 2
     *
     *   SF^2 = pi/2 * (2/2)^2 - 0.04/2
     *        = pi/2 - 0.02
     */
    {
        const auto& bin = result.bin(0);

        assert(bin.count == 2);
        assert(
            std::abs(
                bin.sf_squared -
                (kPi / 2.0 - 0.02)
            ) < 1e-12
        );
    }

    /*
     * Bin 1 (lag 2): pair (0,2).
     *   sum(delta^2) = 4, sum(|delta|) = 2, sum(noise) = 0.02, n = 1
     *
     *   SF^2 = pi/2 * 2^2 - 0.02
     *        = 2*pi - 0.02
     */
    {
        const auto& bin = result.bin(1);

        assert(bin.count == 1);
        assert(
            std::abs(
                bin.sf_squared -
                (2.0 * kPi - 0.02)
            ) < 1e-12
        );
    }


    /*
     * No-noise estimators.
     */
    const auto no_noise =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrderNoNoise
        );

    // Bin 0: 2 / 2 = 1
    assert(
        std::abs(
            no_noise.bin(0).sf_squared - 1.0
        ) < 1e-12
    );

    // Bin 1: 4 / 1 = 4
    assert(
        std::abs(
            no_noise.bin(1).sf_squared - 4.0
        ) < 1e-12
    );

    const auto mad_no_noise =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::MeanAbsoluteDeviationNoNoise
        );

    // Bin 0: pi/2
    assert(
        std::abs(
            mad_no_noise.bin(0).sf_squared -
            (kPi / 2.0)
        ) < 1e-12
    );

    // Bin 1: 2*pi
    assert(
        std::abs(
            mad_no_noise.bin(1).sf_squared -
            (2.0 * kPi)
        ) < 1e-12
    );
}

} // namespace


int main()
{
    test_simple_linear_signal();
    test_noise_correction();
    test_bin_boundaries();
    test_methods();

    return 0;
}