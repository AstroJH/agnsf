#pragma once

#include <string>
#include <vector>

namespace agnsf {
namespace io {

/**
 * A light-curve path together with an optional source redshift.
 *
 * The redshift converts observed lags to the rest frame:
 * dt_rest = dt_obs / (1 + z). redshift == 0 means no correction.
 */
struct PathEntry {
    std::string path;
    double redshift = 0.0;
};


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


/**
 * Read a list of (path, redshift) entries from a text file.
 *
 * Each line may contain one or two whitespace-separated columns:
 *
 *   # readshift = 0
 *   /data/lc1.csv
 *
 *   # readshift = 0.15
 *   /data/lc2.csv 0.15
 *
 * A missing redshift column means no correction (z = 0). Blank lines
 * and '#' comments are ignored.
 *
 * @throws std::invalid_argument for malformed lines or redshift <= -1.
 * @throws std::runtime_error    if the file cannot be opened.
 */
std::vector<PathEntry> read_path_list_with_redshift(
    const std::string& path
);

} // namespace io
} // namespace agnsf
