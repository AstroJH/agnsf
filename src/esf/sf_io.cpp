#include <esf/sf_io.hpp>

#include <stdexcept>
#include <utility>
#include <vector>

namespace agnsf {
namespace esf {

namespace {

std::vector<agnsf::LightCurve> load_curves(
    const std::vector<std::string>& paths,
    const agnsf::io::ColumnNames& columns
)
{
    std::vector<agnsf::LightCurve> curves;
    curves.reserve(paths.size());

    for (const auto& path : paths) {
        curves.emplace_back(
            agnsf::io::read_light_curve(path, columns)
        );
    }

    return curves;
}

} // namespace


SFResult sf_from_file(
    const std::string& path,
    const LagBins& bins,
    SFMethod method,
    const agnsf::io::ColumnNames& columns,
    const UncertaintyConfig& config
)
{
    const agnsf::LightCurve light_curve =
        agnsf::io::read_light_curve(path, columns);

    SFCalculator calculator;

    return calculator.calculate(
        light_curve,
        bins,
        method,
        config
    );
}


SFResult pooled_sf_from_files(
    const std::vector<std::string>& paths,
    const LagBins& bins,
    SFMethod method,
    const agnsf::io::ColumnNames& columns,
    const UncertaintyConfig& config
)
{
    PooledESFCalculator calculator;

    return calculator.calculate(
        load_curves(paths, columns),
        bins,
        method,
        config
    );
}


SFResult pooled_sf_from_path_list(
    const std::string& path_list_file,
    const LagBins& bins,
    SFMethod method,
    const agnsf::io::ColumnNames& columns,
    const UncertaintyConfig& config
)
{
    return pooled_sf_from_files(
        agnsf::io::read_path_list(path_list_file),
        bins,
        method,
        columns,
        config
    );
}


SFResult ensemble_sf_from_files(
    const std::vector<std::string>& paths,
    const LagBins& bins,
    SFMethod sf_method,
    SFEnsembleCalculator::Method method,
    const agnsf::io::ColumnNames& columns,
    const UncertaintyConfig& config
)
{
    SFEnsembleCalculator calculator;

    return calculator.calculate(
        load_curves(paths, columns),
        bins,
        sf_method,
        method,
        config
    );
}


SFResult ensemble_sf_from_path_list(
    const std::string& path_list_file,
    const LagBins& bins,
    SFMethod sf_method,
    SFEnsembleCalculator::Method method,
    const agnsf::io::ColumnNames& columns,
    const UncertaintyConfig& config
)
{
    return ensemble_sf_from_files(
        agnsf::io::read_path_list(path_list_file),
        bins,
        sf_method,
        method,
        columns,
        config
    );
}


void write_sf_result(
    const std::string& path,
    const LagBins& bins,
    const SFResult& result
)
{
    if (result.size() != bins.size()) {
        throw std::invalid_argument(
            "SFResult size does not match the number of lag bins"
        );
    }

    std::vector<double> lag;
    std::vector<double> count;
    std::vector<double> sf_squared;
    std::vector<double> sf;

    lag.reserve(bins.size());
    count.reserve(bins.size());
    sf_squared.reserve(bins.size());
    sf.reserve(bins.size());

    const auto& edges = bins.edges();

    for (std::size_t i = 0; i < bins.size(); ++i) {
        lag.push_back(
            (edges[i] + edges[i + 1]) / 2.0
        );

        count.push_back(
            static_cast<double>(result.bin(i).count)
        );

        sf_squared.push_back(result.bin(i).sf_squared);
        sf.push_back(result.bin(i).sf);
    }

    agnsf::io::write_table(
        path,
        {"lag", "count", "sf_squared", "sf"},
        {lag, count, sf_squared, sf}
    );
}

} // namespace esf
} // namespace agnsf
