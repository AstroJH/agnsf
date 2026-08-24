#include <cassert>
#include <cmath>
#include <vector>

#include <esf/lag_bins.hpp>
#include <esf/light_curve.hpp>
#include <esf/pooled_esf_calculator.hpp>

namespace {

void test_pooled_pairs()
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

    esf::PooledESFCalculator calculator;

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

} // namespace


int main()
{
    test_pooled_pairs();

    return 0;
}