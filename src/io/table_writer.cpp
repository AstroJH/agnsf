#include <io/table_writer.hpp>

#include <fstream>
#include <iomanip>
#include <stdexcept>

namespace agnsf {
namespace io {

void write_table(
    const std::string& path,
    const std::vector<std::string>& headers,
    const std::vector<std::vector<double>>& columns
)
{
    if (headers.size() != columns.size()) {
        throw std::invalid_argument(
            "number of headers must match the number of columns"
        );
    }

    std::size_t nrows = 0;

    if (!columns.empty()) {
        nrows = columns[0].size();
    }

    for (const auto& column : columns) {
        if (column.size() != nrows) {
            throw std::invalid_argument(
                "all columns must have the same length"
            );
        }
    }

    std::ofstream output(path);

    if (!output) {
        throw std::runtime_error(
            "cannot open output file '" + path + "'"
        );
    }

    output << "#";

    for (const auto& header : headers) {
        output << " " << header;
    }

    output << "\n";

    output << std::setprecision(17);

    for (std::size_t i = 0; i < nrows; ++i) {
        for (std::size_t j = 0; j < columns.size(); ++j) {
            if (j != 0) {
                output << " ";
            }

            output << columns[j][i];
        }

        output << "\n";
    }
}

} // namespace io
} // namespace agnsf
