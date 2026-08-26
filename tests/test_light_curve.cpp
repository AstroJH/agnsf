#include <cassert>
#include <cstdint>
#include <stdexcept>

#include <core/light_curve.hpp>

namespace {

void test_basic()
{
    agnsf::LightCurve data(
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
        agnsf::LightCurve data(
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


void test_view_from_light_curve()
{
    agnsf::LightCurve data(
        {0.0, 1.0, 2.0},
        {10.0, 11.0, 12.0},
        {0.1, 0.1, 0.2}
    );

    const agnsf::LightCurveView view =
        data.view();

    assert(view.size() == data.size());

    // A view must point at the owning LightCurve's buffers.
    assert(view.time_data() == data.time_data());
    assert(view.value_data() == data.value_data());
    assert(view.error_data() == data.error_data());

    assert(view.time_address() ==
        reinterpret_cast<std::uintptr_t>(data.time_data()));
    assert(view.value_address() ==
        reinterpret_cast<std::uintptr_t>(data.value_data()));
    assert(view.error_address() ==
        reinterpret_cast<std::uintptr_t>(data.error_data()));

    // The view reflects the underlying data.
    assert(view.time_data()[1] == 1.0);
    assert(view.value_data()[2] == 12.0);
    assert(view.error_data()[0] == 0.1);
}


void test_view_from_raw_arrays()
{
    const double time[] = {0.0, 1.0, 2.0};
    const double value[] = {5.0, 6.0, 7.0};
    const double error[] = {0.01, 0.01, 0.02};

    const agnsf::LightCurveView view(
        time,
        value,
        error,
        3
    );

    assert(view.size() == 3);
    assert(view.time_data() == time);
    assert(view.value_data() == value);
    assert(view.error_data() == error);
}

} // namespace


int main()
{
    test_basic();
    test_size_mismatch();
    test_view_from_light_curve();
    test_view_from_raw_arrays();

    return 0;
}
