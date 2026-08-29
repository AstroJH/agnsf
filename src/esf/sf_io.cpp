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


/**
 * Light curves together with their per-source redshifts (metadata;
 * the redshifts are passed to the calculators, never stored in the
 * LightCurve objects).
 */
struct LoadedCurves {
    std::vector<agnsf::LightCurve> curves;
    std::vector<double> redshifts;
};


/**
 * Load observed-frame light curves and keep their per-source
 * redshifts. No rest-frame time arrays are created here.
 */
LoadedCurves load_curves_with_redshifts(
    const std::vector<agnsf::io::PathEntry>& entries,
    const agnsf::io::ColumnNames& columns
)
{
    LoadedCurves result;

    result.curves.reserve(entries.size());
    result.redshifts.reserve(entries.size());

    for (const auto& entry : entries) {
        result.curves.emplace_back(
            agnsf::io::read_light_curve(entry.path, columns)
        );

        result.redshifts.push_back(entry.redshift);
    }

    return result;
}

} // namespace


SFResult sf_from_file(
    const std::string& path,
    const LagBins& bins,
    SFMethod method,
    const agnsf::io::ColumnNames& columns,
    const UncertaintyConfig& config,
    double redshift
)
{
    const agnsf::LightCurve light_curve =
        agnsf::io::read_light_curve(path, columns);

    SFCalculator calculator;

    return calculator.calculate(
        light_curve,
        bins,
        method,
        config,
        redshift
    );
}


SFResult pooled_sf_from_files(
    const std::vector<std::string>& paths,
    const LagBins& bins,
    SFMethod method,
    const agnsf::io::ColumnNames& columns,
    const UncertaintyConfig& config,
    double redshift
)
{
    PooledESFCalculator calculator;

    return calculator.calculate(
        load_curves(paths, columns),
        bins,
        method,
        config,
        redshift
    );
}


SFResult pooled_sf_from_files(
    const std::vector<std::string>& paths,
    const std::vector<double>& redshifts,
    const LagBins& bins,
    SFMethod method,
    const agnsf::io::ColumnNames& columns,
    const UncertaintyConfig& config
)
{
    if (redshifts.size() != paths.size()) {
        throw std::invalid_argument(
            "redshifts size must match the number of light curves"
        );
    }

    PooledESFCalculator calculator;

    return calculator.calculate(
        load_curves(paths, columns),
        bins,
        method,
        config,
        redshifts
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
    const LoadedCurves loaded =
        load_curves_with_redshifts(
            agnsf::io::read_path_list_with_redshift(path_list_file),
            columns
        );

    PooledESFCalculator calculator;

    return calculator.calculate(
        loaded.curves,
        bins,
        method,
        config,
        loaded.redshifts
    );
}


SFResult ensemble_sf_from_files(
    const std::vector<std::string>& paths,
    const LagBins& bins,
    SFMethod sf_method,
    SFEnsembleCalculator::Method method,
    const agnsf::io::ColumnNames& columns,
    const UncertaintyConfig& config,
    double redshift
)
{
    SFEnsembleCalculator calculator;

    return calculator.calculate(
        load_curves(paths, columns),
        bins,
        sf_method,
        method,
        config,
        redshift
    );
}


SFResult ensemble_sf_from_files(
    const std::vector<std::string>& paths,
    const std::vector<double>& redshifts,
    const LagBins& bins,
    SFMethod sf_method,
    SFEnsembleCalculator::Method method,
    const agnsf::io::ColumnNames& columns,
    const UncertaintyConfig& config
)
{
    if (redshifts.size() != paths.size()) {
        throw std::invalid_argument(
            "redshifts size must match the number of light curves"
        );
    }

    SFEnsembleCalculator calculator;

    return calculator.calculate(
        load_curves(paths, columns),
        bins,
        sf_method,
        method,
        config,
        redshifts
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
    const LoadedCurves loaded =
        load_curves_with_redshifts(
            agnsf::io::read_path_list_with_redshift(path_list_file),
            columns
        );

    SFEnsembleCalculator calculator;

    return calculator.calculate(
        loaded.curves,
        bins,
        sf_method,
        method,
        config,
        loaded.redshifts
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
