#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <io/table_writer.hpp>

namespace fs = std::filesystem;

int main()
{
    const fs::path dir =
        fs::temp_directory_path() / "agnsf_test_table_writer";

    fs::remove_all(dir);
    fs::create_directories(dir);

    const fs::path file = dir / "table.txt";

    agnsf::io::write_table(
        file.string(),
        {"lag", "count", "sf_squared", "sf"},
        {
            {0.5, 1.5, 2.5},
            {3.0, 5.0, 1.0},
            {0.98, 3.98, 1.2345678901234567},
            {0.989949, 1.994993, 1.1111111111111112}
        }
    );

    std::ifstream input(file);
    assert(input);

    std::string line;

    assert(std::getline(input, line));
    assert(line == "# lag count sf_squared sf");

    std::vector<std::vector<double>> rows;

    while (std::getline(input, line)) {
        std::istringstream iss(line);
        std::vector<double> row;
        double value;

        while (iss >> value) {
            row.push_back(value);
        }

        assert(row.size() == 4);
        rows.push_back(row);
    }

    assert(rows.size() == 3);

    assert(std::abs(rows[0][0] - 0.5) < 1e-15);
    assert(std::abs(rows[1][0] - 1.5) < 1e-15);
    assert(std::abs(rows[2][0] - 2.5) < 1e-15);

    assert(std::abs(rows[2][2] - 1.2345678901234567) < 1e-15);

    // Mismatched column lengths must be rejected.
    bool thrown = false;

    try {
        agnsf::io::write_table(
            (dir / "bad.txt").string(),
            {"a", "b"},
            {{1.0, 2.0}, {1.0}}
        );
    }
    catch (const std::invalid_argument&) {
        thrown = true;
    }

    assert(thrown);

    fs::remove_all(dir);

    return 0;
}
