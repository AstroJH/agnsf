#include "bindings/common.hpp"

PYBIND11_MODULE(_agnsf, m)
{
    m.doc() = "Astronomical structure function analysis";

    bind_core(m);
    bind_uncertainty(m);
    bind_sf(m);
    bind_esf(m);
    bind_io(m);
    bind_timedelay(m);
}
