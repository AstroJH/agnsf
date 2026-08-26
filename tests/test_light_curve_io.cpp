#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#include <fitsio.h>

#include <io/light_curve.hpp>

namespace fs = std::filesystem;

namespace {

const fs::path kTempDir =
    fs::temp_directory_path() / "agnsf_test_light_curve_io";


void write_file(
    const fs::path& path,
    const std::string& content
)
{
    std::ofstream out(path);
    assert(out);
    out << content;
}


void check_close(
    double actual,
    double expected
)
{
    assert(
        std::abs(actual - expected) < 1e-12
    );
}


void verify_data(
    const agnsf::LightCurve& data,
    std::size_t expected_size
)
{
    assert(data.size() == expected_size);
    assert(data.time().size() == expected_size);
    assert(data.value().size() == expected_size);
    assert(data.error().size() == expected_size);

    for (std::size_t i = 0; i < expected_size; ++i) {
        check_close(data.time()[i], static_cast<double>(i));
        check_close(data.value()[i], static_cast<double>(i) + 100.0);
        check_close(data.error()[i], 0.1);
    }
}


void test_csv_with_header()
{
    const fs::path path = kTempDir / "with_header.csv";

    write_file(
        path,
        "time,value,error\n"
        "0.0,100.0,0.1\n"
        "1.0,101.0,0.1\n"
        "2.0,102.0,0.1\n"
        "3.0,103.0,0.1\n"
    );

    const agnsf::LightCurve data =
        agnsf::io::read_light_curve(path.string());

    verify_data(data, 4);
}


void test_csv_without_header()
{
    const fs::path path = kTempDir / "no_header.csv";

    write_file(
        path,
        "0.0,100.0,0.1\n"
        "1.0,101.0,0.1\n"
        "2.0,102.0,0.1\n"
    );

    const agnsf::LightCurve data =
        agnsf::io::read_light_curve(path.string());

    verify_data(data, 3);
}


void test_csv_custom_columns()
{
    const fs::path path = kTempDir / "custom_columns.csv";

    write_file(
        path,
        "mjd,mag,emag,flag\n"
        "0.0,100.0,0.1,0\n"
        "1.0,101.0,0.1,1\n"
        "2.0,102.0,0.1,0\n"
    );

    agnsf::io::ColumnNames columns;
    columns.time = "mjd";
    columns.value = "mag";
    columns.error = "emag";

    const agnsf::LightCurve data =
        agnsf::io::read_light_curve(path.string(), columns);

    verify_data(data, 3);
}


void test_csv_missing_column()
{
    const fs::path path = kTempDir / "missing_column.csv";

    write_file(
        path,
        "time,value\n"
        "0.0,100.0\n"
        "1.0,101.0\n"
    );

    bool thrown = false;

    try {
        agnsf::io::read_light_curve(path.string());
    }
    catch (const std::invalid_argument&) {
        thrown = true;
    }

    assert(thrown);
}


void test_unsupported_format()
{
    bool thrown = false;

    try {
        agnsf::io::read_light_curve(
            (kTempDir / "data.txt").string()
        );
    }
    catch (const std::runtime_error&) {
        thrown = true;
    }

    assert(thrown);
}


void write_test_fits(
    const fs::path& path,
    const char* const column_names[3]
)
{
    fitsfile* fptr = nullptr;
    int status = 0;

    fits_create_file(
        &fptr,
        ("!" + path.string()).c_str(),
        &status
    );
    assert(status == 0);

    fits_create_tbl(
        fptr,
        BINARY_TBL,
        0,
        3,
        const_cast<char**>(column_names),
        const_cast<char**>(
            (const char* const[]){"D", "D", "D"}
        ),
        const_cast<char**>(
            (const char* const[]){"", "", ""}
        ),
        "TEST",
        &status
    );
    assert(status == 0);

    double time[] = {0.0, 1.0, 2.0, 3.0};
    double value[] = {100.0, 101.0, 102.0, 103.0};
    double error[] = {0.1, 0.1, 0.1, 0.1};

    fits_write_col(
        fptr, TDOUBLE, 1, 1, 1, 4, time, &status
    );
    assert(status == 0);

    fits_write_col(
        fptr, TDOUBLE, 2, 1, 1, 4, value, &status
    );
    assert(status == 0);

    fits_write_col(
        fptr, TDOUBLE, 3, 1, 1, 4, error, &status
    );
    assert(status == 0);

    fits_close_file(fptr, &status);
    assert(status == 0);
}


void test_fits()
{
    const fs::path path = kTempDir / "curve.fits";

    const char* names[] = {"TIME", "VALUE", "ERROR"};
    write_test_fits(path, names);

    const agnsf::LightCurve data =
        agnsf::io::read_light_curve(path.string());

    verify_data(data, 4);
}


void test_fits_custom_columns()
{
    const fs::path path = kTempDir / "curve_custom.fits";

    const char* names[] = {"MJD", "MAG", "EMAG"};
    write_test_fits(path, names);

    agnsf::io::ColumnNames columns;
    columns.time = "mjd";
    columns.value = "mag";
    columns.error = "emag";

    const agnsf::LightCurve data =
        agnsf::io::read_light_curve(path.string(), columns);

    verify_data(data, 4);
}


void test_fits_missing_column()
{
    const fs::path path = kTempDir / "curve_missing.fits";

    const char* names[] = {"TIME", "VALUE"};
    fitsfile* fptr = nullptr;
    int status = 0;

    fits_create_file(
        &fptr,
        ("!" + path.string()).c_str(),
        &status
    );
    assert(status == 0);

    fits_create_tbl(
        fptr,
        BINARY_TBL,
        0,
        2,
        const_cast<char**>(names),
        const_cast<char**>(
            (const char* const[]){"D", "D"}
        ),
        const_cast<char**>(
            (const char* const[]){"", ""}
        ),
        "TEST",
        &status
    );
    assert(status == 0);

    double time[] = {0.0, 1.0};
    double value[] = {100.0, 101.0};

    fits_write_col(fptr, TDOUBLE, 1, 1, 1, 2, time, &status);
    assert(status == 0);

    fits_write_col(fptr, TDOUBLE, 2, 1, 1, 2, value, &status);
    assert(status == 0);

    fits_close_file(fptr, &status);
    assert(status == 0);

    bool thrown = false;

    try {
        agnsf::io::read_light_curve(path.string());
    }
    catch (const std::invalid_argument&) {
        thrown = true;
    }

    assert(thrown);
}

} // namespace


int main()
{
    fs::remove_all(kTempDir);
    fs::create_directories(kTempDir);

    test_csv_with_header();
    test_csv_without_header();
    test_csv_custom_columns();
    test_csv_missing_column();
    test_unsupported_format();
    test_fits();
    test_fits_custom_columns();
    test_fits_missing_column();

    fs::remove_all(kTempDir);

    return 0;
}
