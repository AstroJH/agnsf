#include <io/light_curve.hpp>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <fitsio.h>

namespace agnsf {
namespace io {

namespace {

std::string to_lower(std::string text)
{
    std::transform(
        text.begin(),
        text.end(),
        text.begin(),
        [](unsigned char c) {
            return static_cast<char>(std::tolower(c));
        }
    );

    return text;
}


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


bool ends_with(
    const std::string& text,
    const std::string& suffix
)
{
    return (
        text.size() >= suffix.size() &&
        text.compare(
            text.size() - suffix.size(),
            suffix.size(),
            suffix
        ) == 0
    );
}


bool try_parse_double(
    const std::string& text,
    double& value
)
{
    if (text.empty()) {
        return false;
    }

    try {
        std::size_t consumed = 0;
        value = std::stod(text, &consumed);
        return consumed == text.size();
    }
    catch (...) {
        return false;
    }
}


/**
 * Split one CSV line into trimmed fields.
 *
 * Handles common CSV quoting rules, including quoted fields,
 * commas inside quoted fields, and escaped double quotes (`""`).
 *
 * This function is intended for parsing a single line and does not
 * perform strict CSV validation.
 *
 * Examples (CSV input -> parsed fields):
 *
 *   Simple input:
 *     1,2,3
 *       -> {"1", "2", "3"}
 *
 *   Quoted field containing a comma:
 *     1,"hello,world",3
 *       -> {"1", "hello,world", "3"}
 *
 *   Escaped double quotes:
 *     "say ""hello"""
 *       -> {"say \"hello\""}
 *
 *   Invalid CSV that is still parsed:
 *     1,"abc,3
 *       -> {"1", "abc,3"}
 *
 *     1,"abc"def,3
 *       -> {"1", "abcdef", "3"}
 *
 *   Valid CSV not supported by this single-line parser:
 *     1,"hello
 *     world",3
 *       -> Not supported (fields spanning multiple lines)
 */
std::vector<std::string> split_csv(
    const std::string& line
)
{
    std::vector<std::string> fields;
    std::string field;

    bool in_quotes = false; // Is inside a quoted field?

    for (std::size_t i = 0; i < line.size(); ++i) {
        const char c = line[i];

        if (in_quotes) {
            if (c == '"') {
                // In CSV, a pair of double quotes inside a quoted field
                // represents a literal double quote.
                if (i + 1 < line.size() && line[i + 1] == '"') {
                    field.push_back('"');
                    ++i; // Skip the second quote.
                } else {

                    // A single quote closes the quoted field.
                    in_quotes = false;
                }
            } else {

                // Characters inside a quoted field are kept as-is,
                // including commas.
                field.push_back(c);
            }
        } else if (c == '"') {
            in_quotes = true; // Start of a quoted field.
        } else if (c == ',') {
            // A comma outside quotes separates two fields.
            fields.push_back(trim(field));
            field.clear();
        } else {
            field.push_back(c);
        }
    }

    // Add the final field, which is not followed by a comma.
    fields.push_back(trim(field));
    return fields;
}


std::size_t find_column(
    const std::vector<std::string>& names,
    const std::string& wanted
)
{
    const std::string target = to_lower(trim(wanted));

    for (std::size_t i = 0; i < names.size(); ++i) {
        if (to_lower(trim(names[i])) == target) {
            return i;
        }
    }

    throw std::invalid_argument(
        "column '" + wanted + "' not found"
    );
}


std::vector<double> parse_fields(
    const std::vector<std::string>& fields,
    const std::vector<std::size_t>& indices
)
{
    std::vector<double> values;
    values.reserve(indices.size());

    for (const std::size_t index : indices) {
        if (index >= fields.size()) {
            throw std::invalid_argument(
                "row has too few columns"
            );
        }

        double value = 0.0;

        if (!try_parse_double(fields[index], value)) {
            throw std::invalid_argument(
                "cannot parse numeric value '" +
                fields[index] + "'"
            );
        }

        values.push_back(value);
    }

    return values;
}


agnsf::LightCurve read_csv(
    const std::string& path,
    const ColumnNames& columns
)
{
    std::ifstream input(path);
    if (!input) {
        throw std::runtime_error(
            "cannot open CSV file '" + path + "'"
        );
    }

    std::string line;
    if (!std::getline(input, line)) {
        throw std::invalid_argument(
            "CSV file '" + path + "' is empty"
        );
    }

    const std::vector<std::string> first = split_csv(line);

    bool has_header = false;

    for (const auto& field : first) {
        double dummy = 0.0;

        if (!try_parse_double(field, dummy)) {
            has_header = true;
            break;
        }
    }

    std::vector<std::size_t> indices(3);

    if (has_header) {
        indices[0] = find_column(first, columns.time);
        indices[1] = find_column(first, columns.value);
        indices[2] = find_column(first, columns.error);
    } else {
        // No header: use the first three columns in the requested order.
        indices = {0, 1, 2};
    }

    std::vector<double> time_values;
    std::vector<double> value_values;
    std::vector<double> error_values;

    if (!has_header) {
        const std::vector<double> values =
            parse_fields(first, indices);

        time_values.push_back(values[0]);
        value_values.push_back(values[1]);
        error_values.push_back(values[2]);
    }

    while (std::getline(input, line)) {
        if (trim(line).empty()) {
            continue;
        }

        const std::vector<double> values =
            parse_fields(split_csv(line), indices);

        time_values.push_back(values[0]);
        value_values.push_back(values[1]);
        error_values.push_back(values[2]);
    }

    if (time_values.empty()) {
        throw std::invalid_argument(
            "CSV file '" + path + "' contains no data rows"
        );
    }

    return agnsf::LightCurve(
        std::move(time_values),
        std::move(value_values),
        std::move(error_values)
    );
}


void check_fits_status(
    int status,
    const std::string& context
)
{
    if (status == 0) {
        return;
    }

    char message[FLEN_ERRMSG] = {0};
    fits_get_errstatus(status, message);

    throw std::runtime_error(
        context + ": " + message
    );
}


agnsf::LightCurve read_fits(
    const std::string& path,
    const ColumnNames& columns
)
{
    fitsfile* fptr = nullptr;
    int status = 0;

    fits_open_file(
        &fptr,
        path.c_str(),
        READONLY,
        &status
    );

    check_fits_status(
        status,
        "cannot open FITS file '" + path + "'"
    );

    // Locate the first table HDU (binary or ASCII).
    int hdutype = 0;
    status = 0;
    fits_get_hdu_type(fptr, &hdutype, &status);
    check_fits_status(status, "reading FITS HDU type");

    while (hdutype != BINARY_TBL && hdutype != ASCII_TBL) {
        status = 0;
        fits_movrel_hdu(fptr, 1, &hdutype, &status);

        if (status != 0) {
            fits_close_file(fptr, &status);
            throw std::invalid_argument(
                "no table HDU found in FITS file '" + path + "'"
            );
        }
    }

    // Get the dimensions of the table.
    int ncols = 0;
    long nrows = 0;
    status = 0;
    fits_get_num_cols(fptr, &ncols, &status);
    fits_get_num_rows(fptr, &nrows, &status);
    check_fits_status(status, "reading FITS table dimensions");

    if (nrows <= 0) {
        fits_close_file(fptr, &status);
        throw std::invalid_argument(
            "FITS table '" + path + "' contains no rows"
        );
    }

    // Read FITS table column names from the TTYPE keywords.
    std::vector<std::string> names(
        static_cast<std::size_t>(ncols)
    );

    for (int i = 1; i <= ncols; ++i) {
        status = 0;

        char keyword[FLEN_KEYWORD];
        std::snprintf(keyword, sizeof(keyword), "TTYPE%d", i);

        char value[FLEN_VALUE] = {0};
        fits_read_key(
            fptr,
            TSTRING,
            keyword,
            value,
            nullptr,
            &status
        );

        if (status == KEY_NO_EXIST) {
            // Keep an empty name for columns without a TTYPE keyword.
            status = 0;
            names[static_cast<std::size_t>(i - 1)] = "";
        } else {
            check_fits_status(status, "reading FITS column names");
            names[static_cast<std::size_t>(i - 1)] = value;
        }
    }

    // Find a column by name, ignoring case and surrounding whitespace.
    const auto column_index =
        [&names](const std::string& wanted) -> int {
            const std::string target = to_lower(trim(wanted));

            for (int i = 0; i < static_cast<int>(names.size()); ++i) {
                if (to_lower(trim(names[static_cast<std::size_t>(i)])) == target) {
                    return i + 1; // CFITSIO column indices are 1-based.
                }
            }

            throw std::invalid_argument(
                "FITS column '" + wanted + "' not found"
            );
        };

    const int time_column = column_index(columns.time);
    const int value_column = column_index(columns.value);
    const int error_column = column_index(columns.error);

    // Read a numeric FITS column into a vector of doubles.
    // FITS null values are represented as NaN.
    const auto read_column =
        [fptr, nrows](int column) -> std::vector<double> {
            std::vector<double> out(
                static_cast<std::size_t>(nrows)
            );

            int status = 0;
            int anynul = 0;
            double nulval =
                std::numeric_limits<double>::quiet_NaN();

            fits_read_col(
                fptr,
                TDOUBLE,
                column,
                1,
                1,
                nrows,
                &nulval,
                out.data(),
                &anynul,
                &status
            );

            check_fits_status(status, "reading FITS column");

            return out;
        };

    // Construct the LightCurve from the selected columns.
    agnsf::LightCurve data(
        read_column(time_column),
        read_column(value_column),
        read_column(error_column)
    );

    status = 0;
    fits_close_file(fptr, &status);

    return data;
}

} // namespace


agnsf::LightCurve read_light_curve(
    const std::string& path,
    const ColumnNames& columns
)
{
    const std::string lower = to_lower(path);

    if (ends_with(lower, ".csv")) {
        return read_csv(path, columns);
    }

    if (
        ends_with(lower, ".fits") ||
        ends_with(lower, ".fit") ||
        ends_with(lower, ".fits.gz") ||
        ends_with(lower, ".fit.gz")
    ) {
        return read_fits(path, columns);
    }

    throw std::runtime_error(
        "unsupported light-curve file format for '" +
        path + "' (supported: .csv, .fits)"
    );
}

} // namespace io
} // namespace agnsf
