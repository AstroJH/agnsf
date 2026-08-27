#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <core/light_curve.hpp>
#include <esf/pooled_esf_calculator.hpp>
#include <esf/sf_calculator.hpp>
#include <esf/sf_ensemble_calculator.hpp>
#include <esf/sf_io.hpp>
#include <io/path_list.hpp>

namespace fs = std::filesystem;

namespace {

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


void check_close(double actual, double expected)
{
    assert(close(actual, expected));
}


agnsf::LightCurve make_curve()
{
    return agnsf::LightCurve(
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 1.0, 1.0, 2.0},
        {0.2, 0.2, 0.2, 0.2}
    );
}


agnsf::LightCurve to_rest_frame(
    const agnsf::LightCurve& curve,
    double redshift
)
{
    std::vector<double> time = curve.time();
    const double factor = 1.0 + redshift;

    for (double& value : time) {
        value /= factor;
    }

    return agnsf::LightCurve(
        time,
        curve.value(),
        curve.error()
    );
}


void check_results_equal(
    const agnsf::esf::SFResult& actual,
    const agnsf::esf::SFResult& expected
)
{
    assert(actual.size() == expected.size());

    for (std::size_t i = 0; i < expected.size(); ++i) {
        assert(actual.bin(i).count == expected.bin(i).count);

        if (expected.bin(i).count == 0) {
            assert(std::isnan(actual.bin(i).sf_squared));
            assert(std::isnan(actual.bin(i).sf));
            continue;
        }

        check_close(
            actual.bin(i).sf_squared,
            expected.bin(i).sf_squared
        );
        check_close(actual.bin(i).sf, expected.bin(i).sf);
        // std::cerr
        // << "bin " << i
        // << ": actual count=" << actual.bin(i).count
        // << ", expected count=" << expected.bin(i).count
        // << ", actual sf2=" << actual.bin(i).sf_squared
        // << ", expected sf2=" << expected.bin(i).sf_squared
        // << "\n";
    }
}


void test_sf_redshift_equivalence()
{
    const double z = 0.5;
    const agnsf::esf::LagBins bins({0.0, 0.5, 0.9, 1.5, 2.2});

    agnsf::esf::SFCalculator calculator;

    const agnsf::esf::SFResult with_redshift =
        calculator.calculate(
            make_curve(),
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            {},
            z
        );

    const agnsf::esf::SFResult manual =
        calculator.calculate(
            to_rest_frame(make_curve(), z),
            bins
        );

    check_results_equal(with_redshift, manual);
}


void test_pooled_redshift_equivalence()
{
    const double z = 0.5;
    const agnsf::esf::LagBins bins({0.0, 0.5, 0.9, 1.5, 2.2});

    agnsf::LightCurve second(
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 2.0, 1.0, 3.0},
        {0.2, 0.2, 0.2, 0.2}
    );

    const std::vector<agnsf::LightCurve> data = {
        make_curve(),
        second
    };

    agnsf::esf::PooledESFCalculator calculator;

    const agnsf::esf::SFResult with_redshift =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            {},
            z
        );

    std::vector<agnsf::LightCurve> manual_data;

    for (const auto& curve : data) {
        manual_data.push_back(to_rest_frame(curve, z));
    }

    const agnsf::esf::SFResult manual =
        calculator.calculate(manual_data, bins);

    check_results_equal(with_redshift, manual);
}


void test_ensemble_redshift_equivalence()
{
    const double z = 0.5;
    const agnsf::esf::LagBins bins({0.0, 0.5, 0.9, 1.5, 2.2});

    agnsf::LightCurve second(
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 2.0, 1.0, 3.0},
        {0.2, 0.2, 0.2, 0.2}
    );

    const std::vector<agnsf::LightCurve> data = {
        make_curve(),
        second
    };

    agnsf::esf::SFEnsembleCalculator calculator;

    const agnsf::esf::SFResult with_redshift =
        calculator.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            agnsf::esf::SFEnsembleCalculator::Method::SqrtMeanSquared,
            {},
            z
        );

    std::vector<agnsf::LightCurve> manual_data;

    for (const auto& curve : data) {
        manual_data.push_back(to_rest_frame(curve, z));
    }

    const agnsf::esf::SFResult manual =
        calculator.calculate(
            manual_data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            agnsf::esf::SFEnsembleCalculator::Method::SqrtMeanSquared
        );

    check_results_equal(with_redshift, manual);
}


void test_invalid_redshift()
{
    agnsf::esf::SFCalculator calculator;
    agnsf::esf::PooledESFCalculator pooled;
    agnsf::esf::SFEnsembleCalculator ensemble;

    const agnsf::esf::LagBins bins({0.0, 1.0});

    bool thrown = false;

    try {
        calculator.calculate(make_curve(), bins, {}, {}, -1.0);
    }
    catch (const std::invalid_argument&) {
        thrown = true;
    }

    assert(thrown);

    thrown = false;

    try {
        pooled.calculate(
            std::vector<agnsf::LightCurve>{make_curve()},
            bins,
            {},
            {},
            -2.0
        );
    }
    catch (const std::invalid_argument&) {
        thrown = true;
    }

    assert(thrown);

    thrown = false;

    try {
        ensemble.calculate(
            std::vector<agnsf::LightCurve>{make_curve()},
            bins,
            {},
            {},
            {},
            -2.0
        );
    }
    catch (const std::invalid_argument&) {
        thrown = true;
    }

    assert(thrown);
}


void test_path_list_with_redshift()
{
    const fs::path dir =
        fs::temp_directory_path() / "agnsf_test_redshift";

    fs::remove_all(dir);
    fs::create_directories(dir);

    const fs::path file = dir / "paths.txt";

    {
        std::ofstream out(file);

        out << "# curves\n"
            << (dir / "lc1.csv").string() << "\n"
            << (dir / "lc2.csv").string() << " 0.5\n";
    }

    const std::vector<agnsf::io::PathEntry> entries =
        agnsf::io::read_path_list_with_redshift(file.string());

    assert(entries.size() == 2);
    assert(entries[0].path == (dir / "lc1.csv").string());
    assert(entries[0].redshift == 0.0);
    assert(entries[1].path == (dir / "lc2.csv").string());
    assert(entries[1].redshift == 0.5);

    // Malformed line: more than two columns.
    {
        const fs::path bad = dir / "bad.txt";

        // Close the writer before reading, so the content is flushed.
        {
            std::ofstream out(bad);
            out << "a.csv 0.5 extra\n";
        }

        bool thrown = false;

        try {
            agnsf::io::read_path_list_with_redshift(bad.string());
        }
        catch (const std::invalid_argument&) {
            thrown = true;
        }

        assert(thrown);
    }

    // redshift <= -1 is rejected.
    {
        const fs::path bad = dir / "bad_z.txt";

        {
            std::ofstream out(bad);
            out << "a.csv -1.0\n";
        }

        bool thrown = false;

        try {
            agnsf::io::read_path_list_with_redshift(bad.string());
        }
        catch (const std::invalid_argument&) {
            thrown = true;
        }

        assert(thrown);
    }

    fs::remove_all(dir);
}


void write_csv(
    const fs::path& path,
    const std::vector<double>& values
)
{
    std::ofstream out(path);

    out << "time,value,error\n";

    for (std::size_t i = 0; i < values.size(); ++i) {
        out << i << "," << values[i] << ",0.2\n";
    }
}


void test_pooled_from_path_list_with_z()
{
    const fs::path dir =
        fs::temp_directory_path() / "agnsf_test_redshift_io";

    fs::remove_all(dir);
    fs::create_directories(dir);

    const fs::path lc1 = dir / "lc1.csv";
    const fs::path lc2 = dir / "lc2.csv";

    write_csv(lc1, {0.0, 1.0, 1.0, 2.0});
    write_csv(lc2, {0.0, 2.0, 1.0, 3.0});

    const fs::path list = dir / "paths.txt";

    {
        std::ofstream out(list);

        out << (lc1.string()) << " 0.5\n"
            << (lc2.string()) << " 0.25\n";
    }

    const agnsf::esf::LagBins bins({0.0, 0.5, 0.9, 1.5, 2.2});

    const agnsf::esf::SFResult from_list =
        agnsf::esf::pooled_sf_from_path_list(
            list.string(),
            bins
        );

    // Manual: scale each curve and pool.
    std::vector<agnsf::LightCurve> curves = {
        to_rest_frame(
            agnsf::io::read_light_curve(lc1.string()),
            0.5
        ),
        to_rest_frame(
            agnsf::io::read_light_curve(lc2.string()),
            0.25
        )
    };

    agnsf::esf::PooledESFCalculator calculator;

    const agnsf::esf::SFResult manual =
        calculator.calculate(curves, bins);

    check_results_equal(from_list, manual);

    fs::remove_all(dir);
}


void test_per_curve_redshift_overload()
{
    const agnsf::esf::LagBins bins({0.0, 0.5, 0.9, 1.5, 2.2});

    agnsf::LightCurve second(
        {0.0, 1.0, 2.0, 3.0},
        {0.0, 2.0, 1.0, 3.0},
        {0.2, 0.2, 0.2, 0.2}
    );

    const std::vector<agnsf::LightCurve> data = {
        make_curve(),
        second
    };

    const std::vector<double> redshifts = {0.5, 0.25};

    // Pooled ESF with per-curve redshifts.
    agnsf::esf::PooledESFCalculator pooled;

    const agnsf::esf::SFResult per_curve =
        pooled.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            {},
            redshifts
        );

    const agnsf::esf::SFResult manual_pooled =
        pooled.calculate(
            {
                to_rest_frame(make_curve(), 0.5),
                to_rest_frame(second, 0.25)
            },
            bins
        );

    check_results_equal(per_curve, manual_pooled);

    // Aggregated ESF with per-curve redshifts.
    agnsf::esf::SFEnsembleCalculator ensemble;

    const agnsf::esf::SFResult per_curve_ensemble =
        ensemble.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            agnsf::esf::SFEnsembleCalculator::Method::SqrtMeanSquared,
            {},
            redshifts
        );

    const agnsf::esf::SFResult manual_ensemble =
        ensemble.calculate(
            {
                to_rest_frame(make_curve(), 0.5),
                to_rest_frame(second, 0.25)
            },
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            agnsf::esf::SFEnsembleCalculator::Method::SqrtMeanSquared
        );

    check_results_equal(per_curve_ensemble, manual_ensemble);

    // The per-curve vector must match the number of curves.
    bool thrown = false;

    try {
        pooled.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            {},
            std::vector<double>{0.5}
        );
    }
    catch (const std::invalid_argument&) {
        thrown = true;
    }

    assert(thrown);

    thrown = false;

    try {
        ensemble.calculate(
            data,
            bins,
            agnsf::esf::SFMethod::SecondOrder,
            agnsf::esf::SFEnsembleCalculator::Method::SqrtMeanSquared,
            {},
            std::vector<double>{0.5, 0.25, 0.1}
        );
    }
    catch (const std::invalid_argument&) {
        thrown = true;
    }

    assert(thrown);
}

} // namespace


int main()
{
    test_sf_redshift_equivalence();
    test_pooled_redshift_equivalence();
    test_ensemble_redshift_equivalence();
    test_invalid_redshift();
    test_path_list_with_redshift();
    test_pooled_from_path_list_with_z();
    test_per_curve_redshift_overload();

    return 0;
}
