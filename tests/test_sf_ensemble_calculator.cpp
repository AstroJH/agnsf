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


    return 0;
}