#pragma once

#include <string>
#include <vector>

namespace agnsf {
namespace io {

/**
 * Read a list of file paths from a text file.
 *
 * One path per line. Empty lines and lines whose first non-space
 * character is '#' are ignored.
 *
 * @throws std::runtime_error if the file cannot be opened.
 */
std::vector<std::string> read_path_list(
    const std::string& path
);

} // namespace io
} // namespace agnsf
