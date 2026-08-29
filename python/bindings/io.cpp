#include "common.hpp"
#include <io/light_curve.hpp>
#include <io/path_list.hpp>
#include <io/table_writer.hpp>

void bind_io(py::module_& m)
{
    m.def(
        "read_light_curve",
        [](
            const std::string& path,
            const std::string& time,
            const std::string& value,
            const std::string& error
        )
        {
            agnsf::io::ColumnNames columns;
            columns.time = time;
            columns.value = value;
            columns.error = error;

            return agnsf::io::read_light_curve(
                path,
                columns
            );
        },
        py::arg("path"),
        py::arg("time") = "time",
        py::arg("value") = "value",
        py::arg("error") = "error",
        R"pbdoc(
Read a light curve from a CSV or FITS file.

The format is chosen from the file extension (.csv, .fits, .fit,
.fits.gz, .fit.gz). Columns are located by name for files with a
header; otherwise the first three columns are used.

Returns a LightCurve with time, value, and error.
        )pbdoc"
    );

    m.def(
        "read_path_list",
        [](
            const std::string& path
        )
        {
            return agnsf::io::read_path_list(path);
        },
        py::arg("path"),
        R"pbdoc(
Read a list of file paths from a text file (one per line, '#' comments
and blank lines are ignored).
        )pbdoc"
    );

    m.def(
        "read_path_list_with_redshift",
        [](
            const std::string& path
        )
        {
            const auto entries =
                agnsf::io::read_path_list_with_redshift(path);

            py::list result;

            for (const auto& entry : entries) {
                result.append(
                    py::make_tuple(entry.path, entry.redshift)
                );
            }

            return result;
        },
        py::arg("path"),
        R"pbdoc(
Read a list of (path, redshift) entries from a text file.

Each line may have one or two whitespace-separated columns:

    /data/lc1.csv
    /data/lc2.csv 0.5

A missing redshift column means no correction (z = 0). Returns a list
of (path, redshift) tuples.
        )pbdoc"
    );

    m.def(
        "write_table",
        [](
            const std::string& path,
            const std::vector<std::string>& headers,
            const std::vector<std::vector<double>>& columns
        )
        {
            agnsf::io::write_table(path, headers, columns);
        },
        py::arg("path"),
        py::arg("headers"),
        py::arg("columns"),
        R"pbdoc(
Write columns of numbers to a simple text file.

The first line is '# header_1 header_2 ...' followed by one row per
line with space-separated values.
        )pbdoc"
    );
}
