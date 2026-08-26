#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

#include <esf/sf_calculator.hpp>
#include <esf/sf_ensemble_calculator.hpp>

namespace agnsf {
namespace esf {

namespace {

/**
 * Turn accumulated per-bin statistics into an SFResult.
 *
 * For SqrtMeanSquared, sum holds the sum of per-curve SF^2 values
 * and the result is:
 *
 *   SF^2 = sum / count
 *   SF   = sqrt(SF^2)
 *
 * For MeanSf, sum holds the sum of per-curve SF values and the
 * result is:
 *
 *   SF   = sum / count
 *   SF^2 = SF^2
 *
 * The SF^2 field is kept consistent with SF in both modes so that
 * sf = sqrt(sf_squared) always holds for finite results.
 */
SFResult make_result(
    const std::vector<double>& sum,
    const std::vector<std::size_t>& count,
    SFEnsembleCalculator::Method method
)
{
    std::vector<SFBinResult> results;
    results.reserve(sum.size());

    for (std::size_t i = 0;
         i < sum.size();
         ++i) {

        SFBinResult result;
        result.count = count[i];

        if (count[i] == 0) {

            result.sf_squared =
                std::numeric_limits<double>::quiet_NaN();

            result.sf =
                std::numeric_limits<double>::quiet_NaN();

        } else if (method == SFEnsembleCalculator::Method::MeanSf) {

            result.sf =
                sum[i] /
                static_cast<double>(count[i]);

            result.sf_squared =
                result.sf * result.sf;

        } else {

            result.sf_squared =
                sum[i] /
                static_cast<double>(count[i]);

            result.sf =
                std::sqrt(result.sf_squared);
        }

        results.push_back(result);
    }

    return SFResult(std::move(results));
}


SFResult calculate_ensemble(
    const std::vector<agnsf::LightCurveView>& data,
    const LagBins& bins,
    SFMethod sf_method,
    SFEnsembleCalculator::Method method
)
{
    if (data.empty()) {

        const std::vector<double> sum(
            bins.size(),
            0.0
        );

        const std::vector<std::size_t> count(
            bins.size(),
            0
        );

        return make_result(
            sum,
            count,
            method
        );
    }


    std::size_t n_threads =
        std::thread::hardware_concurrency();

    if (n_threads == 0) {
        n_threads = 1;
    }

    n_threads =
        std::min(n_threads, data.size());


    /*
     * Each thread has its own partial results.
     */
    std::vector<std::vector<double>> thread_sum(
        n_threads,
        std::vector<double>(bins.size(), 0.0)
    );

    std::vector<std::vector<std::size_t>> thread_count(
        n_threads,
        std::vector<std::size_t>(bins.size(), 0)
    );


    const std::size_t chunk_size =
        (data.size() + n_threads - 1) /
        n_threads;


    std::vector<std::thread> threads;
    threads.reserve(n_threads);


    for (std::size_t thread_id = 0;
         thread_id < n_threads;
         ++thread_id) {

        const std::size_t begin =
            thread_id * chunk_size;

        const std::size_t end =
            std::min(
                begin + chunk_size,
                data.size()
            );

        threads.emplace_back(
            [&, thread_id, begin, end, sf_method, method]() {

                SFCalculator sf_calculator;

                for (std::size_t i = begin;
                     i < end;
                     ++i) {

                    const SFResult result =
                        sf_calculator.calculate(
                            data[i],
                            bins,
                            sf_method
                        );

                    for (std::size_t j = 0;
                         j < bins.size();
                         ++j) {

                        /*
                         * SqrtMeanSquared averages per-curve SF^2.
                         * MeanSf averages per-curve SF instead.
                         */
                        const double contribution =
                            method ==
                                SFEnsembleCalculator::Method::MeanSf
                                ? result.bin(j).sf
                                : result.bin(j).sf_squared;

                        if (!std::isfinite(contribution)) {
                            continue;
                        }

                        thread_sum[thread_id][j] +=
                            contribution;

                        ++thread_count[thread_id][j];
                    }
                }
            }
        );
    }


    for (auto& thread : threads) {
        thread.join();
    }


    /*
     * Merge thread-local results.
     */
    std::vector<double> sum(
        bins.size(),
        0.0
    );

    std::vector<std::size_t> count(
        bins.size(),
        0
    );

    for (std::size_t thread_id = 0;
         thread_id < n_threads;
         ++thread_id) {

        for (std::size_t i = 0;
             i < bins.size();
             ++i) {

            sum[i] +=
                thread_sum[thread_id][i];

            count[i] +=
                thread_count[thread_id][i];
        }
    }

    return make_result(
        sum,
        count,
        method
    );
}

} // namespace


SFResult SFEnsembleCalculator::calculate(
    const std::vector<agnsf::LightCurve>& data,
    const LagBins& bins,
    SFMethod sf_method,
    Method method
) const
{
    std::vector<agnsf::LightCurveView> views;
    views.reserve(data.size());

    for (const auto& light_curve : data) {
        views.push_back(light_curve.view());
    }

    return calculate(
        views,
        bins,
        sf_method,
        method
    );
}


SFResult SFEnsembleCalculator::calculate(
    const std::vector<agnsf::LightCurveView>& data,
    const LagBins& bins,
    SFMethod sf_method,
    Method method
) const
{
    return calculate_ensemble(
        data,
        bins,
        sf_method,
        method
    );
}

} // namespace esf
} // namespace agnsf
