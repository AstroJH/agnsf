#include <io/path_list.hpp>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace agnsf {
namespace io {

namespace {

std::string trim(const std::string& text)
{
    auto is_space =
        [](unsigned char c) {
            return std::isspace(c) != 0;
        };

    const auto first =
        std::find_if_not(
            text.begin(),
            text.end(),
            is_space
        );

    const auto last =
        std::find_if_not(
            text.rbegin(),
            text.rend(),
            is_space
        ).base();

    if (first >= last) {
        return std::string();
    }

    return std::string(first, last);
}


std::vector<std::string> split_whitespace(
    const std::string& text
)
{
    std::istringstream stream(text);
    std::vector<std::string> tokens;
    std::string token;

    while (stream >> token) {
        tokens.push_back(token);
    }

    return tokens;
}


double parse_redshift(
    const std::string& token,
    const std::string& line
)
{
    try {
        std::size_t consumed = 0;
        const double redshift = std::stod(token, &consumed);

        if (consumed != token.size()) {
            throw std::invalid_argument(
                "invalid redshift in path-list line '" + line + "'"
            );
        }

        if (redshift <= -1.0) {
            throw std::invalid_argument(
                "redshift must be > -1 in path-list line '" + line + "'"
            );
        }

        return redshift;
    }
    catch (const std::invalid_argument&) {
        throw;
    }
    catch (const std::exception&) {
        throw std::invalid_argument(
            "invalid redshift in path-list line '" + line + "'"
        );
    }
}

} // namespace


std::vector<std::string> read_path_list(
    const std::string& path
)
{
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "cannot open path list file '" + path + "'"
        );
    }

    std::vector<std::string> paths;
    std::string line;

    while (std::getline(input, line)) {
        const std::string trimmed = trim(line);

        if (trimmed.empty() || trimmed[0] == '#') {
            continue;
        }

        // Keep the full trimmed line so paths may contain spaces.
        paths.push_back(trimmed);
    }

    return paths;
}


std::vector<PathEntry> read_path_list_with_redshift(
    const std::string& path
)
{
    std::ifstream input(path);

    if (!input) {
        throw std::runtime_error(
            "cannot open path list file '" + path + "'"
        );
    }

    std::vector<PathEntry> entries;
    std::string line;

    while (std::getline(input, line)) {
        const std::string trimmed = trim(line);

        if (trimmed.empty() || trimmed[0] == '#') {
            continue;
        }

        const std::vector<std::string> tokens =
            split_whitespace(trimmed);

        if (tokens.empty()) {
            continue;
        }

        if (tokens.size() > 2) {
            throw std::invalid_argument(
                "path-list line must have at most two columns "
                "(path [redshift]): '" + line + "'"
            );
        }

        PathEntry entry;
        entry.path = tokens[0];

        if (tokens.size() == 2) {
            entry.redshift = parse_redshift(tokens[1], line);
        }

        entries.push_back(entry);
    }

    return entries;
}

} // namespace io
} // namespace agnsf
