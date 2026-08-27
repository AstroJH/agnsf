// ---------------------------------------------------------------------------
// agnsf — command-line interface for computing SF/ESF from light-curve
// files.
//
// Examples:
//
//   # SF of one CSV file, written to result.txt
//   agnsf --input lc.csv --output result.txt
//
//   # pooled ESF from a path-list file (@ prefix)
//   agnsf --input @list.txt --output esf.txt --esf pooled
//
//   # aggregated ESF with bootstrap sampling uncertainty
//   agnsf -i @list.txt -o esf.txt --esf aggregated \
//         --sampling bootstrap --n-bootstrap 500 --bootstrap-seed 42
//
// Run `agnsf --help` for the full option list.
// ---------------------------------------------------------------------------

#include <algorithm>
#include <cctype>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <unistd.h>  // isatty / STDIN_FILENO

#include <core/light_curve.hpp>
#include <esf/pooled_esf_calculator.hpp>
#include <esf/sf_calculator.hpp>
#include <esf/sf_ensemble_calculator.hpp>
#include <esf/sf_method.hpp>
#include <esf/sf_result.hpp>
#include <esf/sf_uncertainty.hpp>
#include <io/light_curve.hpp>
#include <io/path_list.hpp>
#include <io/table_writer.hpp>

namespace fs = std::filesystem;

namespace {

using agnsf::esf::LagBins;
using agnsf::esf::SFEnsembleCalculator;
using agnsf::esf::SFMethod;
using agnsf::esf::SFResult;
using agnsf::esf::UncertaintyConfig;
using agnsf::esf::UncertaintyMethod;

// --------------------------------------------------------------------------
// Small helpers
// --------------------------------------------------------------------------

std::string trim(const std::string& text)
{
    auto is_space =
        [](unsigned char c) {
            return std::isspace(c) != 0;
        };

    const auto first =
        std::find_if_not(text.begin(), text.end(), is_space);

    const auto last =
        std::find_if_not(text.rbegin(), text.rend(), is_space).base();

    if (first >= last) {
        return std::string();
    }

    return std::string(first, last);
}


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


std::vector<std::string> split(const std::string& text, char delimiter)
{
    std::vector<std::string> parts;
    std::string part;

    for (const char c : text) {
        if (c == delimiter) {
            parts.push_back(trim(part));
            part.clear();
        } else {
            part.push_back(c);
        }
    }

    parts.push_back(trim(part));
    return parts;
}


std::string format_double(double value)
{
    std::ostringstream out;
    out << value;
    return out.str();
}


// --------------------------------------------------------------------------
// Options
// --------------------------------------------------------------------------

enum class EsfKind {
    Pooled,
    Aggregated
};


struct Options {
    bool help = false;

    // Required (prompted when missing).
    std::string input;
    std::string output;

    // Science settings.
    SFMethod method = SFMethod::SecondOrder;
    EsfKind esf = EsfKind::Pooled;
    SFEnsembleCalculator::Method ensemble_method =
        SFEnsembleCalculator::Method::SqrtMeanSquared;
    UncertaintyConfig uncertainty;
    std::string bins_spec;

    // Column names.
    agnsf::io::ColumnNames columns;

    // Output policy.
    bool clobber = false;
};


void print_help(std::ostream& out)
{
    out <<
        "usage: agnsf [options]\n"
        "\n"
        "Compute the structure function (SF) of one light-curve file, or the\n"
        "ensemble structure function (ESF) of a list of light-curve files.\n"
        "\n"
        "Input:\n"
        "  --input <path>, -i <path>   light-curve file (.csv / .fits), or a\n"
        "                             path-list file when prefixed with '@'\n"
        "                             (e.g. --input @list.txt). An '@' input\n"
        "                             computes the ESF.\n"
        "  --output <path>, -o <path>  output text file. Refuses to overwrite\n"
        "                             an existing file unless --clobber.\n"
        "\n"
        "Science settings:\n"
        "  --method <name>             SF estimator:\n"
        "                               second_order (default)\n"
        "                               second_order_no_noise\n"
        "                               mean_absolute_deviation\n"
        "                               mean_absolute_deviation_no_noise\n"
        "  --esf <kind>                ESF type for '@' inputs:\n"
        "                               pooled (default) | aggregated\n"
        "  --ensemble-method <name>    aggregation for aggregated ESF:\n"
        "                               sqrt_mean_squared (default) | mean_sf\n"
        "  --bins <spec>               lag bins:\n"
        "                               linear:min,max,step (default)\n"
        "                               log:min,max,step\n"
        "                               edges:v1,v2,...\n"
        "                             When omitted, a linear grid from 0 to\n"
        "                             the maximum observed lag (20 bins) is\n"
        "                             used, or an interactive prompt is shown.\n"
        "\n"
        "Uncertainty:\n"
        "  --measurement <mode>        off (default) | analytic\n"
        "  --sampling <mode>           off (default) | analytic | jackknife\n"
        "                               | bootstrap\n"
        "  --n-bootstrap <n>           bootstrap replicates (default 100)\n"
        "  --bootstrap-seed <seed>     bootstrap RNG seed (default 0)\n"
        "\n"
        "Column names (CSV / FITS):\n"
        "  --columns <t,v,e>           comma-separated names, default time,\n"
        "                               value,error\n"
        "  --time-column <name>        override the time column name\n"
        "  --value-column <name>       override the value column name\n"
        "  --error-column <name>       override the error column name\n"
        "\n"
        "Output policy:\n"
        "  --clobber                   overwrite an existing output file\n"
        "  -h, --help                  show this help and exit\n"
        "\n"
        "Output format:\n"
        "  # lag count sf_squared sf measurement_lower measurement_upper\n"
        "  #   sampling_lower sampling_upper\n"
        "\n"
        "Missing required options (--input / --output / --bins) are requested\n"
        "interactively instead of failing immediately.\n"
        "\n"
        "Author: Jiahua Wu (ORCID 0009-0003-1518-6186)\n"
        "Contact: jiahua@stu.xmu.edu.cn (temporary)\n"
        "Report issues: https://github.com/AstroJH/agnsf/issues\n";
}


// --------------------------------------------------------------------------
// Argument parsing
// --------------------------------------------------------------------------

[[noreturn]] void fail(const std::string& message)
{
    throw std::runtime_error(message);
}


double parse_double(const std::string& text, const std::string& option)
{
    try {
        std::size_t consumed = 0;
        const double value = std::stod(text, &consumed);

        if (consumed != text.size()) {
            fail("invalid number '" + text + "' for " + option);
        }

        return value;
    }
    catch (const std::invalid_argument&) {
        fail("invalid number '" + text + "' for " + option);
    }
    catch (const std::out_of_range&) {
        fail("number out of range '" + text + "' for " + option);
    }
}


unsigned long parse_unsigned(const std::string& text, const std::string& option)
{
    try {
        std::size_t consumed = 0;
        const unsigned long value = std::stoul(text, &consumed);

        if (consumed != text.size()) {
            fail("invalid integer '" + text + "' for " + option);
        }

        return value;
    }
    catch (const std::exception&) {
        fail("invalid integer '" + text + "' for " + option);
    }
}


SFMethod parse_sf_method(const std::string& text)
{
    const std::string name = to_lower(text);

    if (name == "second_order") return SFMethod::SecondOrder;
    if (name == "second_order_no_noise") return SFMethod::SecondOrderNoNoise;
    if (name == "mean_absolute_deviation") return SFMethod::MeanAbsoluteDeviation;
    if (name == "mean_absolute_deviation_no_noise") return SFMethod::MeanAbsoluteDeviationNoNoise;

    fail(
        "unknown SF method '" + text + "' (see --help)"
    );
}


EsfKind parse_esf_kind(const std::string& text)
{
    const std::string name = to_lower(text);

    if (name == "pooled") return EsfKind::Pooled;
    if (name == "aggregated") return EsfKind::Aggregated;

    fail(
        "unknown ESF kind '" + text + "' (pooled | aggregated)"
    );
}


SFEnsembleCalculator::Method parse_ensemble_method(const std::string& text)
{
    const std::string name = to_lower(text);

    if (name == "sqrt_mean_squared") {
        return SFEnsembleCalculator::Method::SqrtMeanSquared;
    }

    if (name == "mean_sf") {
        return SFEnsembleCalculator::Method::MeanSf;
    }

    fail(
        "unknown ensemble method '" + text +
        "' (sqrt_mean_squared | mean_sf)"
    );
}


UncertaintyMethod parse_uncertainty_method(
    const std::string& text,
    const std::string& option
)
{
    const std::string name = to_lower(text);

    if (name == "off") return UncertaintyMethod::Off;
    if (name == "analytic") return UncertaintyMethod::Analytic;
    if (name == "jackknife") return UncertaintyMethod::Jackknife;
    if (name == "bootstrap") return UncertaintyMethod::Bootstrap;

    fail(
        "unknown " + option + " '" + text +
        "' (off | analytic | jackknife | bootstrap)"
    );
}


Options parse_args(int argc, char* argv[])
{
    Options options;

    const auto need_value =
        [&options, argc, argv](int& i, const std::string& option,
                               const std::string& inline_value) -> std::string {
            if (!inline_value.empty()) {
                return inline_value;
            }

            if (i + 1 >= argc) {
                fail("missing value for " + option);
            }

            ++i;
            return argv[i];
        };

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];

        std::string name = arg;
        std::string inline_value;

        const std::size_t eq = arg.find('=');

        if (arg.rfind("--", 0) == 0 && eq != std::string::npos) {
            name = arg.substr(0, eq);
            inline_value = arg.substr(eq + 1);
        }

        if (name == "-h" || name == "--help") {
            options.help = true;
        } else if (name == "-i" || name == "--input") {
            options.input = need_value(i, name, inline_value);
        } else if (name == "-o" || name == "--output") {
            options.output = need_value(i, name, inline_value);
        } else if (name == "--method") {
            options.method =
                parse_sf_method(need_value(i, name, inline_value));
        } else if (name == "--esf") {
            options.esf =
                parse_esf_kind(need_value(i, name, inline_value));
        } else if (name == "--ensemble-method") {
            options.ensemble_method =
                parse_ensemble_method(need_value(i, name, inline_value));
        } else if (name == "--measurement") {
            options.uncertainty.measurement =
                parse_uncertainty_method(
                    need_value(i, name, inline_value), "--measurement"
                );
        } else if (name == "--sampling") {
            options.uncertainty.sampling =
                parse_uncertainty_method(
                    need_value(i, name, inline_value), "--sampling"
                );
        } else if (name == "--n-bootstrap") {
            options.uncertainty.n_bootstrap =
                parse_unsigned(need_value(i, name, inline_value),
                               "--n-bootstrap");
        } else if (name == "--bootstrap-seed") {
            options.uncertainty.bootstrap_seed =
                static_cast<std::uint32_t>(
                    parse_unsigned(need_value(i, name, inline_value),
                                   "--bootstrap-seed")
                );
        } else if (name == "--bins") {
            options.bins_spec = need_value(i, name, inline_value);
        } else if (name == "--columns") {
            const std::vector<std::string> names =
                split(need_value(i, name, inline_value), ',');

            if (names.size() != 3) {
                fail("--columns expects exactly three names (time,value,error)");
            }

            options.columns.time = names[0];
            options.columns.value = names[1];
            options.columns.error = names[2];
        } else if (name == "--time-column") {
            options.columns.time = need_value(i, name, inline_value);
        } else if (name == "--value-column") {
            options.columns.value = need_value(i, name, inline_value);
        } else if (name == "--error-column") {
            options.columns.error = need_value(i, name, inline_value);
        } else if (name == "--clobber") {
            options.clobber = true;
        } else {
            fail("unknown option '" + arg + "' (see --help)");
        }
    }

    return options;
}


// --------------------------------------------------------------------------
// Interactive prompting
// --------------------------------------------------------------------------

bool stdin_is_interactive()
{
    return isatty(STDIN_FILENO) != 0;
}


std::string prompt_required(const std::string& question)
{
    if (!stdin_is_interactive()) {
        fail(
            "missing required parameter (" + question +
            ") and stdin is not interactive"
        );
    }

    std::cout << question << ": " << std::flush;

    std::string line;

    if (!std::getline(std::cin, line)) {
        fail("no input provided for " + question);
    }

    const std::string trimmed = trim(line);

    if (trimmed.empty()) {
        fail("empty input provided for " + question);
    }

    return trimmed;
}


std::string prompt_with_default(
    const std::string& question,
    const std::string& default_value
)
{
    if (!stdin_is_interactive()) {
        return default_value;
    }

    std::cout << question
              << " [" << default_value << "]: "
              << std::flush;

    std::string line;

    if (!std::getline(std::cin, line)) {
        return default_value;
    }

    const std::string trimmed = trim(line);

    return trimmed.empty() ? default_value : trimmed;
}


// --------------------------------------------------------------------------
// Bins
// --------------------------------------------------------------------------

LagBins parse_bins_spec(const std::string& spec)
{
    // Forms: "linear:min,max,step" | "log:min,max,step" |
    //        "edges:v1,v2,..." | "min,max,step" (linear shorthand).
    std::string kind = "linear";
    std::string body = spec;

    const std::size_t colon = spec.find(':');

    if (colon != std::string::npos) {
        kind = to_lower(spec.substr(0, colon));
        body = spec.substr(colon + 1);
    }

    const std::vector<std::string> parts = split(body, ',');

    if (kind == "linear") {
        if (parts.size() != 3) {
            fail("linear bins need min,max,step");
        }

        return LagBins::linear(
            parse_double(parts[0], "--bins"),
            parse_double(parts[1], "--bins"),
            parse_double(parts[2], "--bins")
        );
    }

    if (kind == "log") {
        if (parts.size() != 3) {
            fail("log bins need min,max,step");
        }

        return LagBins::logarithmic(
            parse_double(parts[0], "--bins"),
            parse_double(parts[1], "--bins"),
            parse_double(parts[2], "--bins")
        );
    }

    if (kind == "edges") {
        if (parts.size() < 2) {
            fail("edges bins need at least two edges");
        }

        std::vector<double> edges;

        for (const auto& part : parts) {
            edges.push_back(parse_double(part, "--bins"));
        }

        return LagBins(edges);
    }

    fail("unknown bins kind '" + kind + "' (linear | log | edges)");
}


LagBins make_bins(
    const Options& options,
    double max_lag,
    bool interactive
)
{
    if (!options.bins_spec.empty()) {
        return parse_bins_spec(options.bins_spec);
    }

    if (max_lag <= 0.0) {
        fail("cannot derive default lag bins: no positive lag in the data");
    }

    // Default: linear grid from 0 to the maximum observed lag.
    const double step = max_lag / 20.0;

    const std::string suggestion =
        "linear:0," + format_double(max_lag) + "," + format_double(step);

    const std::string spec =
        prompt_with_default("lag bins", suggestion);

    return parse_bins_spec(spec);
}


// --------------------------------------------------------------------------
// Input loading
// --------------------------------------------------------------------------

struct InputData {
    bool is_esf = false;
    std::vector<agnsf::LightCurve> curves;
    double max_lag = 0.0;
};


agnsf::LightCurve load_curve(
    const std::string& path,
    const agnsf::io::ColumnNames& columns
)
{
    return agnsf::io::read_light_curve(path, columns);
}


InputData load_input(const Options& options)
{
    InputData data;

    const std::string& input = options.input;

    if (!input.empty() && input[0] == '@') {

        // Path-list file -> ESF.
        data.is_esf = true;

        const std::string path_list_file = input.substr(1);

        const std::vector<std::string> paths =
            agnsf::io::read_path_list(path_list_file);

        if (paths.empty()) {
            fail("path-list file '" + path_list_file + "' contains no paths");
        }

        for (const auto& path : paths) {
            data.curves.push_back(load_curve(path, options.columns));
        }

    } else {

        // Single light-curve file -> SF.
        data.is_esf = false;
        data.curves.push_back(load_curve(input, options.columns));
    }

    // Maximum observed lag, used for the default bin grid.
    for (const auto& curve : data.curves) {
        if (curve.size() < 2) {
            continue;
        }

        const double curve_lag =
            curve.time().back() - curve.time().front();

        data.max_lag = std::max(data.max_lag, curve_lag);
    }

    return data;
}


// --------------------------------------------------------------------------
// Computation
// --------------------------------------------------------------------------

SFResult compute(
    const Options& options,
    const InputData& data,
    const LagBins& bins
)
{
    if (data.is_esf) {

        if (options.esf == EsfKind::Pooled) {

            agnsf::esf::PooledESFCalculator calculator;

            return calculator.calculate(
                data.curves,
                bins,
                options.method,
                options.uncertainty
            );
        }

        agnsf::esf::SFEnsembleCalculator calculator;

        return calculator.calculate(
            data.curves,
            bins,
            options.method,
            options.ensemble_method,
            options.uncertainty
        );
    }

    agnsf::esf::SFCalculator calculator;

    return calculator.calculate(
        data.curves.front(),
        bins,
        options.method,
        options.uncertainty
    );
}


// --------------------------------------------------------------------------
// Settings printing (before execution)
// --------------------------------------------------------------------------

std::string uncertainty_method_name(UncertaintyMethod method)
{
    switch (method) {
        case UncertaintyMethod::Off: return "off";
        case UncertaintyMethod::Analytic: return "analytic";
        case UncertaintyMethod::Jackknife: return "jackknife";
        case UncertaintyMethod::Bootstrap: return "bootstrap";
    }

    return "?";
}


std::string sf_method_name(SFMethod method)
{
    switch (method) {
        case SFMethod::SecondOrder: return "second_order";
        case SFMethod::SecondOrderNoNoise: return "second_order_no_noise";
        case SFMethod::MeanAbsoluteDeviation: return "mean_absolute_deviation";
        case SFMethod::MeanAbsoluteDeviationNoNoise:
            return "mean_absolute_deviation_no_noise";
    }

    return "?";
}


void print_settings(
    const Options& options,
    const InputData& data,
    const LagBins& bins
)
{
    std::cout << "# agnsf settings\n"
              << "#   input:  " << options.input << "\n"
              << "#   output: " << options.output << "\n";

    if (data.is_esf) {
        std::cout
            << "#   mode:   ESF ("
            << (options.esf == EsfKind::Pooled ? "pooled" : "aggregated")
            << ")\n";
    } else {
        std::cout << "#   mode:   SF\n";
    }

    std::cout
        << "#   method: " << sf_method_name(options.method) << "\n";

    if (data.is_esf && options.esf == EsfKind::Aggregated) {
        std::cout
            << "#   ensemble method: "
            << (options.ensemble_method ==
                        SFEnsembleCalculator::Method::SqrtMeanSquared
                    ? "sqrt_mean_squared"
                    : "mean_sf")
            << "\n";
    }

    std::cout
        << "#   bins:   " << bins.size() << " bins in ["
        << format_double(bins.min()) << ", "
        << format_double(bins.max()) << ")\n"
        << "#   measurement uncertainty: "
        << uncertainty_method_name(options.uncertainty.measurement) << "\n"
        << "#   sampling uncertainty: "
        << uncertainty_method_name(options.uncertainty.sampling) << "\n";

    if (options.uncertainty.sampling == UncertaintyMethod::Bootstrap) {
        std::cout
            << "#   n_bootstrap: " << options.uncertainty.n_bootstrap << "\n"
            << "#   bootstrap_seed: "
            << options.uncertainty.bootstrap_seed << "\n";
    }

    std::cout
        << "#   columns: time='" << options.columns.time
        << "' value='" << options.columns.value
        << "' error='" << options.columns.error << "'\n"
        << "#   clobber: " << (options.clobber ? "yes" : "no") << "\n";
}


// --------------------------------------------------------------------------
// Output
// --------------------------------------------------------------------------

void write_output(
    const std::string& path,
    const LagBins& bins,
    const SFResult& result,
    bool clobber
)
{
    if (fs::exists(path) && !clobber) {
        fail(
            "output file '" + path + "' already exists; "
            "use --clobber to overwrite"
        );
    }

    const auto& edges = bins.edges();

    std::vector<double> lag;
    std::vector<double> count;
    std::vector<double> sf_squared;
    std::vector<double> sf;
    std::vector<double> measurement_lower;
    std::vector<double> measurement_upper;
    std::vector<double> sampling_lower;
    std::vector<double> sampling_upper;

    lag.reserve(bins.size());
    count.reserve(bins.size());
    sf_squared.reserve(bins.size());
    sf.reserve(bins.size());
    measurement_lower.reserve(bins.size());
    measurement_upper.reserve(bins.size());
    sampling_lower.reserve(bins.size());
    sampling_upper.reserve(bins.size());

    for (std::size_t i = 0; i < bins.size(); ++i) {

        const auto& bin = result.bin(i);

        lag.push_back((edges[i] + edges[i + 1]) / 2.0);
        count.push_back(static_cast<double>(bin.count));
        sf_squared.push_back(bin.sf_squared);
        sf.push_back(bin.sf);
        measurement_lower.push_back(bin.measurement.lower);
        measurement_upper.push_back(bin.measurement.upper);
        sampling_lower.push_back(bin.sampling.lower);
        sampling_upper.push_back(bin.sampling.upper);
    }

    agnsf::io::write_table(
        path,
        {
            "lag",
            "count",
            "sf_squared",
            "sf",
            "measurement_lower",
            "measurement_upper",
            "sampling_lower",
            "sampling_upper"
        },
        {
            lag,
            count,
            sf_squared,
            sf,
            measurement_lower,
            measurement_upper,
            sampling_lower,
            sampling_upper
        }
    );
}

} // namespace


int main(int argc, char* argv[])
{
    try {
        Options options = parse_args(argc, argv);

        if (options.help) {
            print_help(std::cout);
            return 0;
        }

        // Required parameters are requested interactively when missing.
        if (options.input.empty()) {
            options.input =
                prompt_required(
                    "input (light-curve file, or @path-list file)"
                );
        }

        if (options.output.empty()) {
            options.output = prompt_required("output file");
        }

        const InputData data = load_input(options);

        const LagBins bins =
            make_bins(options, data.max_lag, stdin_is_interactive());

        // Print the resolved settings before executing.
        print_settings(options, data, bins);

        const SFResult result = compute(options, data, bins);

        write_output(options.output, bins, result, options.clobber);

        std::cout << "# wrote " << options.output << "\n";

        return 0;
    }
    catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << "\n";
        return 1;
    }
}
