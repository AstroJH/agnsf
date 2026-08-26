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


void test_methods()
{
    constexpr double kPi = 3.14159265358979323846;

    esf::BinAccumulator accumulator;

    accumulator.add(1.0, 0.1, 0.1);
    accumulator.add(2.0, 0.1, 0.1);

    /*
     * sum(delta^2) = 1 + 4 = 5
     * sum(|delta|) = 1 + 2 = 3
     * sum(noise)   = 0.02 + 0.02 = 0.04
     * count        = 2
     */
    assert(accumulator.count() == 2);
    assert(std::abs(accumulator.sum_abs_delta() - 3.0) < 1e-12);

    // SecondOrder: (5 - 0.04) / 2 = 2.48
    assert(
        std::abs(
            accumulator.sf_squared(
                esf::SFMethod::SecondOrder
            ) - 2.48
        ) < 1e-12
    );

    // SecondOrderNoNoise: 5 / 2 = 2.5
    assert(
        std::abs(
            accumulator.sf_squared(
                esf::SFMethod::SecondOrderNoNoise
            ) - 2.5
        ) < 1e-12
    );

    // MeanAbsoluteDeviation:
    //   pi/2 * (3/2)^2 - 0.04/2 = pi/2 * 2.25 - 0.02
    const double mad =
        kPi / 2.0 * 2.25 - 0.02;

    assert(
        std::abs(
            accumulator.sf_squared(
                esf::SFMethod::MeanAbsoluteDeviation
            ) - mad
        ) < 1e-12
    );

    assert(
        std::abs(
            accumulator.sf(
                esf::SFMethod::MeanAbsoluteDeviation
            ) - std::sqrt(mad)
        ) < 1e-12
    );

    // MeanAbsoluteDeviationNoNoise:
    //   pi/2 * (3/2)^2 = pi/2 * 2.25
    const double mad_no_noise =
        kPi / 2.0 * 2.25;

    assert(
        std::abs(
            accumulator.sf_squared(
                esf::SFMethod::MeanAbsoluteDeviationNoNoise
            ) - mad_no_noise
        ) < 1e-12
    );

    // Default method must be SecondOrder.
    assert(
        std::abs(
            accumulator.sf_squared() - 2.48
        ) < 1e-12
    );
}

void test_merge()
{
    esf::BinAccumulator left;
    esf::BinAccumulator right;

    left.add(1.0, 0.1, 0.1);
    right.add(-2.0, 0.2, 0.3);

    left.merge(right);

    /*
     * sum(delta^2) = 1 + 4 = 5
     * sum(|delta|) = 1 + 2 = 3
     * sum(noise)   = 0.02 + 0.13 = 0.15
     * count        = 2
     */
    assert(left.count() == 2);
    assert(std::abs(left.sum_delta_squared() - 5.0) < 1e-12);
    assert(std::abs(left.sum_abs_delta() - 3.0) < 1e-12);
    assert(std::abs(left.sum_noise() - 0.15) < 1e-12);

    // SecondOrderNoNoise: 5 / 2 = 2.5
    assert(
        std::abs(
            left.sf_squared(
                esf::SFMethod::SecondOrderNoNoise
            ) - 2.5
        ) < 1e-12
    );
}

} // namespace


int main()
{
    test_noise_correction();
    test_negative_sf_squared();
    test_methods();
    test_merge();

    return 0;
}