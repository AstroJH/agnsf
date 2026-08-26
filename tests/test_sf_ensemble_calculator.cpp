#include <cassert>
#include <cmath>
#include <cstddef>
#include <vector>

#include <esf/lag_bins.hpp>
#include <esf/light_curve.hpp>
#include <esf/light_curve_view.hpp>
#include <esf/sf_calculator.hpp>
#include <esf/sf_ensemble_calculator.hpp>

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

void check_close(
    double actual,
    double expected
)
{
    assert(close(actual, expected));
}

} // namespace



void test_mean_sf_mode()
{
    /*
     * lc_signal:  strong signal -> positive SF everywhere.
     * lc_noise:   pure noise -> negative SF^2, so SF is NaN.
     * lc_short:   only one pair in bin 1 -> bin 0 has no pairs.
     */
    const esf::LightCurve lc_signal(
        {0.0, 1.0, 2.0},
        {0.0, 10.0, 20.0},
        {0.0, 0.0, 0.0}
    );

    const esf::LightCurve lc_noise(
        {0.0, 1.0, 2.0},
        {0.0, 0.0, 0.0},
        {1.0, 1.0, 1.0}
    );

    const esf::LightCurve lc_short(
        {0.0, 2.0},
        {0.0, 4.0},
        {0.0, 0.0}
    );

    const std::vector<esf::LightCurve> data = {
        lc_signal,
        lc_noise,
        lc_short
    };

    const esf::LagBins bins({
        0.0,
        1.5,
        3.0
    });

    esf::SFEnsembleCalculator ensemble_calculator;
    esf::SFCalculator sf_calculator;

    const std::vector<esf::SFResult> individual = {
        sf_calculator.calculate(lc_signal, bins),
        sf_calculator.calculate(lc_noise, bins),
        sf_calculator.calculate(lc_short, bins)
    };

    /*
     * Sanity check the individual SFs used below:
     *
     *   lc_signal: bin0 SF^2=100  SF=10
     *              bin1 SF^2=400  SF=20
     *
     *   lc_noise:  bin0 SF^2=-2   SF=NaN
     *              bin1 SF^2=-2   SF=NaN
     *
     *   lc_short:  bin0 no pairs  -> SF^2=NaN  SF=NaN
     *              bin1 SF^2=16   SF=4
     */
    assert(
        std::abs(
            individual[0].bin(0).sf_squared - 100.0
        ) < 1e-12
    );
    assert(
        std::abs(
            individual[0].bin(1).sf_squared - 400.0
        ) < 1e-12
    );
    assert(
        std::abs(
            individual[1].bin(0).sf_squared + 2.0
        ) < 1e-12
    );
    assert(
        std::isnan(individual[1].bin(0).sf)
    );
    assert(
        std::isnan(individual[2].bin(0).sf_squared)
    );
    assert(
        std::abs(
            individual[2].bin(1).sf_squared - 16.0
        ) < 1e-12
    );


    /*
     * MeanSf mode:
     *
     *   ESF = <SF_k>
     *
     * Only curves with a finite SF contribute.  lc_noise has
     * NaN SF in both bins, so it is excluded from both means.
     */
    const esf::SFResult mean_sf =
        ensemble_calculator.calculate(
            data,
            bins,
            esf::SFMethod::SecondOrder,
            esf::SFEnsembleCalculator::Method::MeanSf
        );

    {
        const auto& bin = mean_sf.bin(0);

        assert(bin.count == 1);
        check_close(bin.sf, 10.0);
        check_close(bin.sf_squared, 100.0);
    }

    {
        const auto& bin = mean_sf.bin(1);

        // (20 + 4) / 2 = 12
        assert(bin.count == 2);
        check_close(bin.sf, 12.0);
        check_close(bin.sf_squared, 144.0);
    }


    /*
     * SqrtMeanSquared mode:
     *
     *   ESF = sqrt(<SF^2>)
     *
     * Negative finite SF^2 values are kept.  lc_short has no
     * pairs in bin 0, so it is excluded there.
     */
    const esf::SFResult sqrt_mean_squared =
        ensemble_calculator.calculate(
            data,
            bins,
            esf::SFMethod::SecondOrder,
            esf::SFEnsembleCalculator::Method::SqrtMeanSquared
        );

    {
        const auto& bin = sqrt_mean_squared.bin(0);

        // (100 + (-2)) / 2 = 49
        assert(bin.count == 2);
        check_close(bin.sf_squared, 49.0);
        check_close(bin.sf, 7.0);
    }

    {
        const auto& bin = sqrt_mean_squared.bin(1);

        // (400 + (-2) + 16) / 3 = 138
        assert(bin.count == 3);
        check_close(bin.sf_squared, 138.0);
        check_close(bin.sf, std::sqrt(138.0));
    }


    /*
     * The default method must be SqrtMeanSquared.
     */
    const esf::SFResult default_result =
        ensemble_calculator.calculate(
            data,
            bins
        );

    for (std::size_t i = 0;
         i < bins.size();
         ++i) {

        assert(
            default_result.bin(i).count ==
            sqrt_mean_squared.bin(i).count
        );

        check_close(
            default_result.bin(i).sf_squared,
            sqrt_mean_squared.bin(i).sf_squared
        );
    }


    /*
     * LightCurveView interface must agree with the owning
     * LightCurve interface for both methods.
     */
    const std::vector<esf::LightCurveView> views = {
        esf::LightCurveView(
            lc_signal.time_data(),
            lc_signal.value_data(),
            lc_signal.error_data(),
            lc_signal.size()
        ),
        esf::LightCurveView(
            lc_noise.time_data(),
            lc_noise.value_data(),
            lc_noise.error_data(),
            lc_noise.size()
        ),
        esf::LightCurveView(
            lc_short.time_data(),
            lc_short.value_data(),
            lc_short.error_data(),
            lc_short.size()
        )
    };

    const esf::SFResult mean_sf_views =
        ensemble_calculator.calculate(
            views,
            bins,
            esf::SFMethod::SecondOrder,
            esf::SFEnsembleCalculator::Method::MeanSf
        );

    for (std::size_t i = 0;
         i < bins.size();
         ++i) {

        assert(
            mean_sf_views.bin(i).count ==
            mean_sf.bin(i).count
        );

        if (mean_sf.bin(i).count == 0) {
            assert(
                std::isnan(mean_sf_views.bin(i).sf)
            );
            continue;
        }

        check_close(
            mean_sf_views.bin(i).sf,
            mean_sf.bin(i).sf
        );

        check_close(
            mean_sf_views.bin(i).sf_squared,
            mean_sf.bin(i).sf_squared
        );
    }


    /*
     * Empty input.
     */
    const std::vector<esf::LightCurve> empty_data;

    const esf::SFResult empty_mean_sf =
        ensemble_calculator.calculate(
            empty_data,
            bins,
            esf::SFMethod::SecondOrder,
            esf::SFEnsembleCalculator::Method::MeanSf
        );

    assert(
        empty_mean_sf.size() == bins.size()
    );

    for (std::size_t i = 0;
         i < empty_mean_sf.size();
         ++i) {

        assert(empty_mean_sf.bin(i).count == 0);
        assert(
            std::isnan(empty_mean_sf.bin(i).sf_squared)
        );
        assert(
            std::isnan(empty_mean_sf.bin(i).sf)
        );
    }
}

void test_sf_method_variants()
{
    constexpr double kPi = 3.14159265358979323846;

    /*
     * Two curves, each with a single pair at lag 1.
     *
     * lc_a: delta = 2  ->  SF^2 = pi/2 * 2^2 = 2*pi
     * lc_b: delta = 4  ->  SF^2 = pi/2 * 4^2 = 8*pi
     */
    const esf::LightCurve lc_a(
        {0.0, 1.0},
        {0.0, 2.0},
        {0.0, 0.0}
    );

    const esf::LightCurve lc_b(
        {0.0, 1.0},
        {0.0, 4.0},
        {0.0, 0.0}
    );

    const std::vector<esf::LightCurve> data = {
        lc_a,
        lc_b
    };

    const esf::LagBins bins({
        0.0,
        1.5
    });

    esf::SFEnsembleCalculator ensemble_calculator;

    /*
     * SqrtMeanSquared with MeanAbsoluteDeviation per-curve SF:
     *
     *   ESF^2 = <SF^2> = (2*pi + 8*pi) / 2 = 5*pi
     */
    const esf::SFResult rms =
        ensemble_calculator.calculate(
            data,
            bins,
            esf::SFMethod::MeanAbsoluteDeviation,
            esf::SFEnsembleCalculator::Method::SqrtMeanSquared
        );

    assert(rms.bin(0).count == 2);
    check_close(
        rms.bin(0).sf_squared,
        5.0 * kPi
    );
    check_close(
        rms.bin(0).sf,
        std::sqrt(5.0 * kPi)
    );

    /*
     * MeanSf with MeanAbsoluteDeviation per-curve SF:
     *
     *   ESF = <SF> = (sqrt(2*pi) + sqrt(8*pi)) / 2
     */
    const esf::SFResult mean =
        ensemble_calculator.calculate(
            data,
            bins,
            esf::SFMethod::MeanAbsoluteDeviation,
            esf::SFEnsembleCalculator::Method::MeanSf
        );

    assert(mean.bin(0).count == 2);
    check_close(
        mean.bin(0).sf,
        (
            std::sqrt(2.0 * kPi) +
            std::sqrt(8.0 * kPi)
        ) / 2.0
    );
    check_close(
        mean.bin(0).sf_squared,
        mean.bin(0).sf * mean.bin(0).sf
    );


    /*
     * SecondOrderNoNoise per-curve SF:
     *
     *   ESF^2 = <SF^2> = (4 + 16) / 2 = 10
     */
    const esf::SFResult no_noise =
        ensemble_calculator.calculate(
            data,
            bins,
            esf::SFMethod::SecondOrderNoNoise,
            esf::SFEnsembleCalculator::Method::SqrtMeanSquared
        );

    assert(no_noise.bin(0).count == 2);
    check_close(no_noise.bin(0).sf_squared, 10.0);
    check_close(no_noise.bin(0).sf, std::sqrt(10.0));
}


int main()
{
    /*
     * Use several light curves with different:
     *
     *   - numbers of observations
     *   - sampling patterns
     *   - variability amplitudes
     *   - photometric errors
     *
     * All times are sorted.
     */

    const esf::LightCurve lc1(
        {
            0.0, 0.8, 1.7, 3.1, 4.0,
            5.6, 7.2, 8.1, 10.0, 12.5
        },
        {
            15.2, 15.5, 15.0, 15.8, 15.3,
            16.1, 15.7, 16.0, 15.4, 16.3
        },
        {
            0.03, 0.04, 0.03, 0.05, 0.04,
            0.05, 0.04, 0.03, 0.05, 0.04
        }
    );

    const esf::LightCurve lc2(
        {
            0.2, 1.0, 2.4, 3.0, 4.8,
            6.0, 7.5, 9.3, 11.0
        },
        {
            14.8, 14.9, 15.4, 15.1, 15.8,
            15.6, 16.2, 15.9, 16.5
        },
        {
            0.05, 0.04, 0.06, 0.05, 0.07,
            0.05, 0.06, 0.05, 0.07
        }
    );

    const esf::LightCurve lc3(
        {
            0.0, 1.5, 2.0, 4.2, 5.0,
            6.8, 8.9, 10.5, 13.0, 15.0,
            17.5
        },
        {
            16.1, 15.7, 16.0, 15.2, 15.5,
            14.9, 15.4, 14.8, 15.1, 14.6,
            15.0
        },
        {
            0.02, 0.03, 0.02, 0.04, 0.03,
            0.04, 0.03, 0.05, 0.04, 0.05,
            0.04
        }
    );

    const esf::LightCurve lc4(
        {
            0.4, 1.1, 2.8, 4.0, 5.7,
            8.0, 10.2, 12.0, 14.5
        },
        {
            15.9, 15.9, 16.0, 15.8, 16.1,
            15.9, 16.2, 16.0, 16.3
        },
        {
            0.06, 0.05, 0.07, 0.06, 0.05,
            0.06, 0.05, 0.07, 0.06
        }
    );

    const std::vector<esf::LightCurve> data = {
        lc1,
        lc2,
        lc3,
        lc4
    };


    /*
     * Several lag bins.
     *
     * The largest bins will naturally contain fewer
     * pairs for some light curves.
     */
    const esf::LagBins bins = esf::LagBins::linear(
        0.5,
        8.5,
        1.0
    );


    /*
     * Calculate the ensemble SF.
     */
    esf::SFEnsembleCalculator ensemble_calculator;

    const esf::SFResult ensemble =
        ensemble_calculator.calculate(
            data,
            bins
        );


    /*
     * Independently calculate the individual SFs.
     */
    esf::SFCalculator sf_calculator;

    const std::vector<esf::SFResult> individual = {
        sf_calculator.calculate(lc1, bins),
        sf_calculator.calculate(lc2, bins),
        sf_calculator.calculate(lc3, bins),
        sf_calculator.calculate(lc4, bins)
    };


    /*
     * Independently construct:
     *
     *     <SF²>
     *
     * and
     *
     *     ESF = sqrt(<SF²>)
     */
    for (std::size_t i = 0;
         i < bins.size();
         ++i) {

        double sum_sf_squared = 0.0;
        std::size_t count = 0;

        for (const auto& result : individual) {

            const double sf_squared =
                result.bin(i).sf_squared;

            if (!std::isfinite(sf_squared)) {
                continue;
            }

            sum_sf_squared += sf_squared;
            ++count;
        }

        assert(
            ensemble.bin(i).count == count
        );

        if (count == 0) {

            assert(
                std::isnan(
                    ensemble.bin(i).sf_squared
                )
            );

            assert(
                std::isnan(
                    ensemble.bin(i).sf
                )
            );

            continue;
        }

        const double expected_sf_squared =
            sum_sf_squared /
            static_cast<double>(count);

        const double expected_sf =
            std::sqrt(expected_sf_squared);

        check_close(
            ensemble.bin(i).sf_squared,
            expected_sf_squared
        );

        check_close(
            ensemble.bin(i).sf,
            expected_sf
        );
    }


    /*
     * Check that the LightCurveView interface produces
     * exactly the same result.
     */
    const std::vector<esf::LightCurveView> views = {
        esf::LightCurveView(
            lc1.time_data(),
            lc1.value_data(),
            lc1.error_data(),
            lc1.size()
        ),
        esf::LightCurveView(
            lc2.time_data(),
            lc2.value_data(),
            lc2.error_data(),
            lc2.size()
        ),
        esf::LightCurveView(
            lc3.time_data(),
            lc3.value_data(),
            lc3.error_data(),
            lc3.size()
        ),
        esf::LightCurveView(
            lc4.time_data(),
            lc4.value_data(),
            lc4.error_data(),
            lc4.size()
        )
    };

    const esf::SFResult ensemble_from_views =
        ensemble_calculator.calculate(
            views,
            bins
        );

    assert(
        ensemble_from_views.size() ==
        ensemble.size()
    );

    for (std::size_t i = 0;
         i < ensemble.size();
         ++i) {

        assert(
            ensemble_from_views.bin(i).count ==
            ensemble.bin(i).count
        );

        if (ensemble.bin(i).count == 0) {
            assert(
                std::isnan(
                    ensemble_from_views.bin(i).sf
                )
            );
            continue;
        }

        check_close(
            ensemble_from_views.bin(i).sf_squared,
            ensemble.bin(i).sf_squared
        );

        check_close(
            ensemble_from_views.bin(i).sf,
            ensemble.bin(i).sf
        );
    }


    /*
     * Empty input.
     */
    const std::vector<esf::LightCurve> empty_data;

    const esf::SFResult empty_result =
        ensemble_calculator.calculate(
            empty_data,
            bins
        );

    assert(
        empty_result.size() == bins.size()
    );

    for (std::size_t i = 0;
         i < empty_result.size();
         ++i) {

        assert(
            empty_result.bin(i).count == 0
        );

        assert(
            std::isnan(
                empty_result.bin(i).sf_squared
            )
        );

        assert(
            std::isnan(
                empty_result.bin(i).sf
            )
        );
    }


    test_mean_sf_mode();
    test_sf_method_variants();


    return 0;
}