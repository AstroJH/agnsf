#pragma once

#include <string>
#include <vector>

namespace agnsf {
namespace io {

/**
 * Write columns of double values to a simple text file.
 *
 * The first line is a header comment:
 *
 *   # header_1 header_2 ...
 *
 * followed by one row per line with space-separated values. Values
 * are written with enough precision to round-trip (17 significant
 * digits).
 *
 * @param headers Column names.
 * @param columns One vector per column; all must share the same length.
 *
 * @throws std::invalid_argument if the column lengths differ.
 * @throws std::runtime_error    if the file cannot be opened.
 */
void write_table(
    const std::string& path,
    const std::vector<std::string>& headers,
    const std::vector<std::vector<double>>& columns
);

} // namespace io
} // namespace agnsf
