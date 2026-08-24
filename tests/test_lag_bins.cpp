#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <esf/lag_bins.hpp>

namespace {

void test_custom_bins()
{
    esf::LagBins bins({
        0.0,
        1.0,
        2.0,
        5.0
    });

    assert(bins.size() == 3);
    assert(bins.min() == 0.0);
    assert(bins.max() == 5.0);

    assert(bins.contains(0.0));
    assert(bins.contains(0.5));
    assert(bins.contains(1.0));
    assert(bins.contains(4.999));

    assert(!bins.contains(5.0));
    assert(!bins.contains(-1.0));
}


void test_index()
{
    esf::LagBins bins({
        0.0,
        1.0,
        2.0,
        5.0
    });

    assert(bins.index(0.0) == 0);
    assert(bins.index(0.999) == 0);

    assert(bins.index(1.0) == 1);
    assert(bins.index(1.999) == 1);

    assert(bins.index(2.0) == 2);
    assert(bins.index(4.999) == 2);
}


void test_try_index()
{
    esf::LagBins bins({
        0.0,
        1.0,
        2.0,
        5.0
    });

    std::size_t index = 999;

    assert(bins.try_index(0.0, index));
    assert(index == 0);

    assert(bins.try_index(1.0, index));
    assert(index == 1);

    assert(bins.try_index(4.999, index));
    assert(index == 2);

    assert(!bins.try_index(5.0, index));
    assert(!bins.try_index(-1.0, index));

    assert(
        !bins.try_index(
            std::numeric_limits<double>::quiet_NaN(),
            index
        )
    );
}


void test_invalid_edges()
{
    bool thrown = false;

    try {
        esf::LagBins bins({
            -1.0,
            1.0
        });
    }
    catch (const std::invalid_argument&) {
        thrown = true;
    }

    assert(thrown);


    thrown = false;

    try {
        esf::LagBins bins({
            0.0,
            1.0,
            1.0
        });
    }
    catch (const std::invalid_argument&) {
        thrown = true;
    }

    assert(thrown);
}


void test_linear()
{
    const auto bins =
        esf::LagBins::linear(
            0.0,
            5.0,
            1.0
        );

    assert(bins.size() == 5);

    assert(bins.edges()[0] == 0.0);
    assert(bins.edges()[1] == 1.0);
    assert(bins.edges()[5] == 5.0);
}


void test_logarithmic()
{
    const auto bins =
        esf::LagBins::logarithmic(
            1.0,
            100.0,
            1.0
        );

    /*
     * log10(1)   = 0
     * log10(10)  = 1
     * log10(100) = 2
     *
     * Therefore:
     *
     * edges = [1, 10, 100]
     * size  = 2
     */

    assert(bins.size() == 2);

    assert(
        std::abs(
            bins.edges()[0] - 1.0
        ) < 1e-12
    );

    assert(
        std::abs(
            bins.edges()[1] - 10.0
        ) < 1e-12
    );

    assert(
        std::abs(
            bins.edges()[2] - 100.0
        ) < 1e-12
    );
}


void test_logarithmic_zero()
{
    bool thrown = false;

    try {
        esf::LagBins::logarithmic(
            0.0,
            100.0,
            10
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
    test_custom_bins();
    test_index();
    test_try_index();
    test_invalid_edges();
    test_linear();
    test_logarithmic();
    test_logarithmic_zero();

    return 0;
}