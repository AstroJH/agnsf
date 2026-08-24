#include <cassert>
#include <cmath>

#include <esf/lag_bins.hpp>
#include <esf/light_curve.hpp>
#include <esf/sf_calculator.hpp>

namespace {

void test_simple_linear_signal()
{
    esf::LightCurve data(
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 0.0, 0.0, 0.0}
    );

    esf::LagBins bins({
        0.0,
        1.5,
        2.5,
        4.0
    });

    esf::SFCalculator calculator;

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
    esf::LightCurve data(
        {0.0, 1.0, 2.0},
        {0.0, 1.0, 2.0},
        {0.1, 0.1, 0.1}
    );

    esf::LagBins bins({
        0.0,
        1.5,
        3.0
    });

    esf::SFCalculator calculator;

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
    esf::LightCurve data(
        {0.0, 1.0, 2.0},
        {0.0, 1.0, 2.0},
        {0.0, 0.0, 0.0}
    );

    esf::LagBins bins({
        0.0,
        1.0,
        2.0,
        3.0
    });

    esf::SFCalculator calculator;

    const auto result =
        calculator.calculate(data, bins);

    // tau = 1 belongs to [1, 2)
    assert(result.bin(0).count == 0);
    assert(result.bin(1).count == 2);

    // tau = 2 belongs to [2, 3)
    assert(result.bin(2).count == 1);
}
} // namespace


int main()
{
    test_simple_linear_signal();
    test_noise_correction();
    test_bin_boundaries();

    return 0;
}