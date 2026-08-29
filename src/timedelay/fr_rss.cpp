#include <timedelay/fr_rss.hpp>

#include <algorithm>
#include <cmath>
#include <random>
#include <vector>

namespace agnsf {
namespace timedelay {

namespace {

/**
 * Build one perturbed / subsampled realization of a light curve.
 *
 * - RSS: sample `n` epochs with replacement, sorted by time;
 * - FR: perturb each sampled flux by N(0, sigma_i^2).
 */
agnsf::LightCurve realization(
    const agnsf::LightCurve& light_curve,
    bool random_subset,
    bool flux_randomization,
    std::mt19937& rng,
    std::normal_distribution<double>& normal
)
{
    const std::size_t n = light_curve.size();

    std::vector<std::size_t> indices(n);

    for (std::size_t i = 0; i < n; ++i) {
        indices[i] = random_subset
            ? static_cast<std::size_t>(rng() % n)
            : i;
    }

    if (random_subset) {
        std::sort(indices.begin(), indices.end());
    }

    std::vector<double> time;
    std::vector<double> value;
    std::vector<double> error;

    time.reserve(n);
    value.reserve(n);
    error.reserve(n);

    for (const std::size_t index : indices) {

        time.push_back(light_curve.time_data()[index]);

        double v = light_curve.value_data()[index];
        const double e = light_curve.error_data()[index];

        if (flux_randomization) {
            v += e * normal(rng);
        }

        value.push_back(v);
        error.push_back(e);
    }

    return agnsf::LightCurve(
        std::move(time),
        std::move(value),
        std::move(error)
    );
}

} // namespace


agnsf::Uncertainty lag_uncertainty(
    const agnsf::LightCurve& continuum,
    const agnsf::LightCurve& response,
    const LagGrid& grid,
    LagEstimate estimate,
    const CrossCorrelationConfig& ccf_config,
    const FRRSSConfig& config
)
{
    agnsf::Uncertainty result;

    if (config.n_realizations == 0) {
        return result;
    }

    std::mt19937 rng(config.seed);
    std::normal_distribution<double> normal(0.0, 1.0);

    std::vector<double> lags;
    lags.reserve(config.n_realizations);

    for (std::size_t r = 0; r < config.n_realizations; ++r) {

        const agnsf::LightCurve continuum_real =
            realization(
                continuum,
                config.random_subset,
                config.flux_randomization,
                rng,
                normal
            );

        const agnsf::LightCurve response_real =
            realization(
                response,
                config.random_subset,
                config.flux_randomization,
                rng,
                normal
            );

        const LagResult cc =
            cross_correlate(
                continuum_real,
                response_real,
                grid,
                ccf_config
            );

        const double lag = lag_value(cc, estimate);

        if (lag == lag) {
            lags.push_back(lag);
        }
    }

    if (lags.size() < 2) {
        return result;
    }

    std::sort(lags.begin(), lags.end());

    const std::size_t lower_index =
        static_cast<std::size_t>(
            0.16 * static_cast<double>(lags.size() - 1)
        );

    const std::size_t upper_index =
        static_cast<std::size_t>(
            0.84 * static_cast<double>(lags.size() - 1)
        );

    result.lower = lags[lower_index];
    result.upper = lags[upper_index];

    return result;
}

} // namespace timedelay
} // namespace agnsf
