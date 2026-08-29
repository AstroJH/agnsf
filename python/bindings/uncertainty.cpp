#include "common.hpp"

void bind_uncertainty(py::module_& m)
{
    py::class_<agnsf::Uncertainty>(m, "Uncertainty")
    .def_readonly(
        "lower",
        &agnsf::Uncertainty::lower
    )
    .def_readonly(
        "upper",
        &agnsf::Uncertainty::upper
    )
    .def_property_readonly(
        "estimated",
        &agnsf::Uncertainty::estimated
    )
    .def(
        "__repr__",
        [](const agnsf::Uncertainty& uncertainty)
        {
            std::ostringstream out;
            out << "Uncertainty(lower="
                << uncertainty.lower
                << ", upper="
                << uncertainty.upper
                << ")";
            return out.str();
        }
    );


    // ------------------------------------------------------------------
    // UncertaintyMethod
    // ------------------------------------------------------------------
    py::enum_<agnsf::UncertaintyMethod>(
        m,
        "UncertaintyMethod"
    )
    .value(
        "Off",
        agnsf::UncertaintyMethod::Off
    )
    .value(
        "Analytic",
        agnsf::UncertaintyMethod::Analytic
    )
    .value(
        "MonteCarlo",
        agnsf::UncertaintyMethod::MonteCarlo
    )
    .value(
        "Jackknife",
        agnsf::UncertaintyMethod::Jackknife
    )
    .value(
        "Bootstrap",
        agnsf::UncertaintyMethod::Bootstrap
    );


    // ------------------------------------------------------------------
    // UncertaintyConfig
    // ------------------------------------------------------------------
    py::class_<agnsf::UncertaintyConfig>(
        m,
        "UncertaintyConfig"
    )
    .def(
        py::init(
            [](
                agnsf::UncertaintyMethod measurement,
                agnsf::UncertaintyMethod within,
                agnsf::UncertaintyMethod sampling,
                std::size_t n_bootstrap,
                std::uint32_t bootstrap_seed
            )
            {
                agnsf::UncertaintyConfig config;

                config.measurement = measurement;
                config.within = within;
                config.sampling = sampling;
                config.n_bootstrap = n_bootstrap;
                config.bootstrap_seed = bootstrap_seed;

                return config;
            }
        ),
        py::arg("measurement") =
            agnsf::UncertaintyMethod::Off,
        py::arg("within") =
            agnsf::UncertaintyMethod::Off,
        py::arg("sampling") =
            agnsf::UncertaintyMethod::Off,
        py::arg("n_bootstrap") =
            std::size_t{100},
        py::arg("bootstrap_seed") =
            std::uint32_t{0}
    )
    .def_readwrite(
        "measurement",
        &agnsf::UncertaintyConfig::measurement
    )
    .def_readwrite(
        "within",
        &agnsf::UncertaintyConfig::within
    )
    .def_readwrite(
        "sampling",
        &agnsf::UncertaintyConfig::sampling
    )
    .def_readwrite(
        "n_bootstrap",
        &agnsf::UncertaintyConfig::n_bootstrap
    )
    .def_readwrite(
        "bootstrap_seed",
        &agnsf::UncertaintyConfig::bootstrap_seed
    )
    .def(
        "__repr__",
        [](const agnsf::UncertaintyConfig& config)
        {
            const auto name =
                [](agnsf::UncertaintyMethod method)
                {
                    switch (method) {
                        case agnsf::UncertaintyMethod::Off:
                            return "Off";
                        case agnsf::UncertaintyMethod::Analytic:
                            return "Analytic";
                        case agnsf::UncertaintyMethod::MonteCarlo:
                            return "MonteCarlo";
                        case agnsf::UncertaintyMethod::Jackknife:
                            return "Jackknife";
                        case agnsf::UncertaintyMethod::Bootstrap:
                            return "Bootstrap";
                    }
                    return "?";
                };

            std::ostringstream out;
            out << "UncertaintyConfig(measurement="
                << name(config.measurement)
                << ", within="
                << name(config.within)
                << ", sampling="
                << name(config.sampling)
                << ", n_bootstrap="
                << config.n_bootstrap
                << ", bootstrap_seed="
                << config.bootstrap_seed
                << ")";
            return out.str();
        }
    );
}
