#pragma once

#include <string>

#include <core/light_curve.hpp>

namespace agnsf {
namespace io {

/**
 * Column names used to locate the time/value/error columns.
 *
 * The default column names:
 *
 *   time  = "time"
 *   value = "value"
 *   error = "error"
 */
struct ColumnNames {
    std::string time = "time";
    std::string value = "value";
    std::string error = "error";
};


/**
 * Read a light curve from a CSV or FITS file.
 *
 * The format is chosen from the file extension:
 *
 *   comma-separated text: .csv
 *   FITS table (binary or ASCII): .fits/.fit/.fits.gz/.fit.gz
 *
 * For CSV files with a header row, columns are located by name (case-insensitive).
 * Without a header, the first three columns are
 * used in the order given by `columns`.
 *
 * For FITS files, the first table extension (or the primary HDU when
 * it is a table) is read, and columns are located by TTYPE name
 * (case-insensitive).
 *
 * @throws std::invalid_argument for missing columns or malformed data.
 * @throws std::runtime_error    for unreadable files or unsupported formats.
 */
agnsf::LightCurve read_light_curve(
    const std::string& path,
    const ColumnNames& columns = {}
);

} // namespace io
} // namespace agnsf
