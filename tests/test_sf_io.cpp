#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include <esf/sf_io.hpp>

namespace fs = std::filesystem;

namespace {

const fs::path kTempDir =
    fs::temp_directory_path() / "agnsf_test_sf_io";

void write_file(
    const fs::path& path,
    const std::string& content
)
{
    std::ofstream out(path);
    assert(out);
    out << content;
}


bool close(
    double a,
    double b,
    double rtol = 1e-12,
    double atol = 1e-12
)
{
    return std::abs(a - b) <=
        atol + rtol * std::abs(b);
}


void check_close(
    double actual,
    double expected
)
{
    assert(close(actual, expected));
}


agnsf::LightCurve make_lc1()
{
    return agnsf::LightCurve(
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 0.0, 0.0, 0.0}
    );
}


agnsf::LightCurve make_lc2()
{
    return agnsf::LightCurve(
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 2.0, 4.0, 6.0},
        {0.0, 0.0, 0.0, 0.0}
    );
}


void test_sf_from_file()
{
    const fs::path csv = kTempDir / "lc1.csv";

    write_file(
        csv,
        "time,value,error\n"
        "0.0,0.0,0.0\n"
        "1.0,1.0,0.0\n"
        "2.0,2.0,0.0\n"
        "3.0,3.0,0.0\n"
    );

    const agnsf::esf::LagBins bins({0.0, 1.5, 2.5, 4.0});

    const agnsf::esf::SFResult from_file =
        agnsf::esf::sf_from_file(csv.string(), bins);

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult direct =
        calculator.calculate(make_lc1(), bins);

    assert(from_file.size() == direct.size());

    for (std::size_t i = 0; i < direct.size(); ++i) {
        assert(from_file.bin(i).count == direct.bin(i).count);
        check_close(
            from_file.bin(i).sf_squared,
            direct.bin(i).sf_squared
        );
        check_close(
            from_file.bin(i).sf,
            direct.bin(i).sf
        );
    }
}


void test_pooled_from_files()
{
    const fs::path csv1 = kTempDir / "lc1.csv";
    const fs::path csv2 = kTempDir / "lc2.csv";

    write_file(
        csv2,
        "time,value,error\n"
        "0.0,0.0,0.0\n"
        "1.0,2.0,0.0\n"
        "2.0,4.0,0.0\n"
        "3.0,6.0,0.0\n"
    );

    const agnsf::esf::LagBins bins({0.0, 1.5, 2.5, 4.0});

    const agnsf::esf::SFResult from_files =
        agnsf::esf::pooled_sf_from_files(
            {csv1.string(), csv2.string()},
            bins
        );

    agnsf::esf::PooledESFCalculator calculator;

    const agnsf::esf::SFResult direct =
        calculator.calculate(
            std::vector<agnsf::LightCurve>{
                make_lc1(),
                make_lc2()
            },
            bins
        );

    assert(from_files.size() == direct.size());

    for (std::size_t i = 0; i < direct.size(); ++i) {
        assert(from_files.bin(i).count == direct.bin(i).count);
        check_close(
            from_files.bin(i).sf_squared,
            direct.bin(i).sf_squared
        );
    }

    // Manual pooled check:
    //   bin0: (1 + 4) / 2 = 2.5
    //   bin1: (4 + 16) / 2 = 10
    //   bin2: (9 + 36) / 2 = 22.5
    check_close(from_files.bin(0).sf_squared, 2.5);
    check_close(from_files.bin(1).sf_squared, 10.0);
    check_close(from_files.bin(2).sf_squared, 22.5);
}


void test_pooled_from_path_list()
{
    const fs::path path_list =
        kTempDir / "paths.txt";

    write_file(
        path_list,
        "# light curves\n"
        + (kTempDir / "lc1.csv").string() + "\n"
        + (kTempDir / "lc2.csv").string() + "\n"
    );

    const agnsf::esf::LagBins bins({0.0, 1.5, 2.5, 4.0});

    const agnsf::esf::SFResult from_list =
        agnsf::esf::pooled_sf_from_path_list(
            path_list.string(),
            bins
        );

    const agnsf::esf::SFResult from_files =
        agnsf::esf::pooled_sf_from_files(
            {
                (kTempDir / "lc1.csv").string(),
                (kTempDir / "lc2.csv").string()
            },
            bins
        );

    assert(from_list.size() == from_files.size());

    for (std::size_t i = 0; i < from_files.size(); ++i) {
        assert(from_list.bin(i).count == from_files.bin(i).count);
        check_close(
            from_list.bin(i).sf_squared,
            from_files.bin(i).sf_squared
        );
    }
}


void test_ensemble_from_files()
{
    const agnsf::esf::LagBins bins({0.0, 1.5, 2.5, 4.0});

    const agnsf::esf::SFResult from_files =
        agnsf::esf::ensemble_sf_from_files(
            {
                (kTempDir / "lc1.csv").string(),
                (kTempDir / "lc2.csv").string()
            },
            bins
        );

    agnsf::esf::SFEnsembleCalculator calculator;

    const agnsf::esf::SFResult direct =
        calculator.calculate(
            std::vector<agnsf::LightCurve>{
                make_lc1(),
                make_lc2()
            },
            bins
        );

    assert(from_files.size() == direct.size());

    for (std::size_t i = 0; i < direct.size(); ++i) {
        assert(from_files.bin(i).count == direct.bin(i).count);
        check_close(
            from_files.bin(i).sf_squared,
            direct.bin(i).sf_squared
        );
        check_close(
            from_files.bin(i).sf,
            direct.bin(i).sf
        );
    }
}


void test_write_sf_result()
{
    const agnsf::esf::LagBins bins({0.0, 1.5, 2.5, 4.0});

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult result =
        calculator.calculate(make_lc1(), bins);

    const fs::path output =
        kTempDir / "result.txt";

    agnsf::esf::write_sf_result(
        output.string(),
        bins,
        result
    );

    std::ifstream input(output);
    assert(input);

    std::string line;

    assert(std::getline(input, line));
    assert(line == "# lag count sf_squared sf");

    // Expected:
    //   lag = 0.75, 2.0, 3.25
    //   count = 3, 2, 1
    //   sf_squared = 1, 4, 9
    //   sf = 1, 2, 3
    const double expected[3][4] = {
        {0.75, 3.0, 1.0, 1.0},
        {2.0, 2.0, 4.0, 2.0},
        {3.25, 1.0, 9.0, 3.0}
    };

    std::size_t row = 0;

    while (std::getline(input, line)) {
        assert(row < 3);

        std::istringstream iss(line);
        double values[4] = {0, 0, 0, 0};

        for (double& value : values) {
            assert(iss >> value);
        }

        for (int j = 0; j < 4; ++j) {
            check_close(
                values[j],
                expected[row][j]
            );
        }

        ++row;
    }

    assert(row == 3);
}

} // namespace


int main()
{
    fs::remove_all(kTempDir);
    fs::create_directories(kTempDir);

    test_sf_from_file();
    test_pooled_from_files();
    test_pooled_from_path_list();
    test_ensemble_from_files();
    test_write_sf_result();

    fs::remove_all(kTempDir);

    return 0;
}
