#include <cassert>
#include <stdexcept>

#include <esf/light_curve.hpp>

namespace {

void test_basic()
{
    esf::LightCurve data(
        {0.0, 1.0, 2.0},
        {10.0, 11.0, 12.0},
        {0.1, 0.1, 0.2}
    );

    assert(data.size() == 3);

    assert(data.time()[0] == 0.0);
    assert(data.value()[1] == 11.0);
    assert(data.error()[2] == 0.2);
}


void test_size_mismatch()
{
    bool thrown = false;

    try {
        esf::LightCurve data(
            {0.0, 1.0},
            {10.0},
            {0.1, 0.1}
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
    test_basic();
    test_size_mismatch();

    return 0;
}