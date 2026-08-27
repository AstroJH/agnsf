#pragma once

#include <string>
#include <vector>

#include <io/light_curve.hpp>
#include <io/path_list.hpp>
#include <io/table_writer.hpp>

#include <esf/lag_bins.hpp>
#include <esf/pooled_esf_calculator.hpp>
#include <esf/sf_calculator.hpp>
#include <esf/sf_ensemble_calculator.hpp>
#include <esf/sf_method.hpp>
#include <esf/sf_result.hpp>
#include <esf/sf_uncertainty.hpp>

namespace agnsf {
namespace esf {

/**
 * Read one light curve from a file and calculate its SF.
 */
SFResult sf_from_file(
    const std::string& path,
    const LagBins& bins,
    SFMethod method = SFMethod::SecondOrder,
    const agnsf::io::ColumnNames& columns = {},
    const UncertaintyConfig& config = {},
    double redshift = 0.0
);


/**
 * Read several light curves from files and calculate the pooled ESF.
 */
SFResult pooled_sf_from_files(
    const std::vector<std::string>& paths,
    const LagBins& bins,
    SFMethod method = SFMethod::SecondOrder,
    const agnsf::io::ColumnNames& columns = {},
    const UncertaintyConfig& config = {},
    double redshift = 0.0
);


/**
 * Read a list of light-curve paths from a text file and calculate
 * the pooled ESF.
 *
 * Each line may carry an optional second column with the source
 * redshift (see agnsf::io::read_path_list_with_redshift); per-curve
 * redshifts are applied before computing.
 */
SFResult pooled_sf_from_path_list(
    const std::string& path_list_file,
    const LagBins& bins,
    SFMethod method = SFMethod::SecondOrder,
    const agnsf::io::ColumnNames& columns = {},
    const UncertaintyConfig& config = {}
);


/**
 * Read several light curves from files and calculate the aggregated
 * ESF.
 */
SFResult ensemble_sf_from_files(
    const std::vector<std::string>& paths,
    const LagBins& bins,
    SFMethod sf_method = SFMethod::SecondOrder,
    SFEnsembleCalculator::Method method =
        SFEnsembleCalculator::Method::SqrtMeanSquared,
    const agnsf::io::ColumnNames& columns = {},
    const UncertaintyConfig& config = {},
    double redshift = 0.0
);


/**
 * Read a list of light-curve paths from a text file and calculate
 * the aggregated ESF.
 *
 * Each line may carry an optional second column with the source
 * redshift; per-curve redshifts are applied before computing.
 */
SFResult ensemble_sf_from_path_list(
    const std::string& path_list_file,
    const LagBins& bins,
    SFMethod sf_method = SFMethod::SecondOrder,
    SFEnsembleCalculator::Method method =
        SFEnsembleCalculator::Method::SqrtMeanSquared,
    const agnsf::io::ColumnNames& columns = {},
    const UncertaintyConfig& config = {}
);


/**
 * Write an SFResult to a simple text file.
 *
 * Columns:
 *
 *   # lag count sf_squared sf
 *
 * The lag column represents the lag associated with each bin.
 *
 * TODO: Consider recording the bin edges explicitly instead of a
 *       single representative lag value.
 */
void write_sf_result(
    const std::string& path,
    const LagBins& bins,
    const SFResult& result
);

} // namespace esf
} // namespace agnsf
