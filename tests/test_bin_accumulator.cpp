#include <cassert>
#include <cmath>

#include <esf/bin_accumulator.hpp>

namespace {

void test_noise_correction()
{
    esf::BinAccumulator accumulator;

    accumulator.add(
        1.0,
        0.1,
        0.1
    );

    accumulator.add(
        2.0,
        0.1,
        0.1
    );

    /*
     * First pair:
     *
     * 1^2 - (0.1^2 + 0.1^2)
     * = 0.98
     *
     * Second pair:
     *
     * 2^2 - 0.02
     * = 3.98
     *
     * Mean:
     *
     * (0.98 + 3.98) / 2 = 2.48
     */
    assert(accumulator.count() == 2);

    assert(
        std::abs(
            accumulator.sf_squared() - 2.48
        ) < 1e-12
    );

    assert(
        std::abs(
            accumulator.sf() - std::sqrt(2.48)
        ) < 1e-12
    );
}

void test_negative_sf_squared()
{
    esf::BinAccumulator accumulator;

    accumulator.add(
        0.1,
        1.0,
        1.0
    );

    assert(accumulator.count() == 1);

    assert(
        std::abs(
            accumulator.sf_squared() + 1.99
        ) < 1e-12
    );

    assert(
        std::isnan(accumulator.sf())
    );
}

} // namespace


int main()
{
    test_noise_correction();
    test_negative_sf_squared();

    return 0;
}