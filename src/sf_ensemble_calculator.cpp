#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <thread>
#include <utility>
#include <vector>

#include <esf/sf_calculator.hpp>
#include <esf/sf_ensemble_calculator.hpp>

namespace esf {

namespace {

LightCurveView make_view(
    const LightCurve& data
)
{
    return LightCurveView(
        data.time_data(),
        data.value_data(),
        data.error_data(),
        data.size()
    );
}


SFResult make_result(
    const std::vector<double>& sum_sf_squared,
    const std::vector<std::size_t>& count
)
{
    std::vector<SFBinResult> results;
    results.reserve(sum_sf_squared.size());

    for (std::size_t i = 0;
         i < sum_sf_squared.size();
         ++i) {

        SFBinResult result;
        result.count = count[i];

        if (count[i] == 0) {

            result.sf_squared =
                std::numeric_limits<double>::quiet_NaN();

            result.sf =
                std::numeric_limits<double>::quiet_NaN();

        } else {

            result.sf_squared =
                sum_sf_squared[i] /
                static_cast<double>(count[i]);

            result.sf =
                std::sqrt(result.sf_squared);
        }

        results.push_back(result);
    }

    return SFResult(std::move(results));
}


SFResult calculate_ensemble(
    const std::vector<LightCurveView>& data,
    const LagBins& bins
)
{
    if (data.empty()) {

        const std::vector<double> sum_sf_squared(
            bins.size(),
            0.0
        );

        const std::vector<std::size_t> count(
            bins.size(),
            0
        );

        return make_result(
            sum_sf_squared,
            count
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
            [&, thread_id, begin, end]() {

                SFCalculator sf_calculator;

                for (std::size_t i = begin;
                     i < end;
                     ++i) {

                    const SFResult result =
                        sf_calculator.calculate(
                            data[i],
                            bins
                        );

                    for (std::size_t j = 0;
                         j < bins.size();
                         ++j) {

                        const double sf_squared =
                            result.bin(j).sf_squared;

                        if (!std::isfinite(sf_squared)) {
                            continue;
                        }

                        thread_sum[thread_id][j] +=
                            sf_squared;

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
    std::vector<double> sum_sf_squared(
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

            sum_sf_squared[i] +=
                thread_sum[thread_id][i];

            count[i] +=
                thread_count[thread_id][i];
        }
    }

    return make_result(
        sum_sf_squared,
        count
    );
}

} // namespace


SFResult SFEnsembleCalculator::calculate(
    const std::vector<LightCurve>& data,
    const LagBins& bins
) const
{
    std::vector<LightCurveView> views;
    views.reserve(data.size());

    for (const auto& light_curve : data) {
        views.push_back(
            make_view(light_curve)
        );
    }

    return calculate(
        views,
        bins
    );
}


SFResult SFEnsembleCalculator::calculate(
    const std::vector<LightCurveView>& data,
    const LagBins& bins
) const
{
    return calculate_ensemble(
        data,
        bins
    );
}

} // namespace esf