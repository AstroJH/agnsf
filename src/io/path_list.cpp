#include <io/path_list.hpp>

#include <algorithm>
#include <cctype>
#include <fstream>
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

        paths.push_back(trimmed);
    }

    return paths;
}

} // namespace io
} // namespace agnsf
