#include <timedelay/cross_correlation.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace agnsf {
namespace timedelay {

namespace {

bool is_finite(double value)
{
    return value == value;
}


std::size_t LagGrid_size(const LagGrid& grid)
{
    if (grid.step <= 0.0 || grid.max <= grid.min) {
        return 0;
    }

    const std::size_t n =
        static_cast<std::size_t>(
            (grid.max - grid.min) / grid.step + 0.5
        );

    return n + 1;
}


std::vector<double> make_taus(const LagGrid& grid)
{
    std::vector<double> taus;
    taus.reserve(LagGrid_size(grid));

    for (std::size_t k = 0; k < LagGrid_size(grid); ++k) {
        taus.push_back(grid.min + static_cast<double>(k) * grid.step);
    }

    return taus;
}


double pearson(
    const std::vector<double>& x,
    const std::vector<double>& y
)
{
    const double n = static_cast<double>(x.size());

    if (n == 0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    double mean_x = 0.0;
    double mean_y = 0.0;

    for (std::size_t i = 0; i < x.size(); ++i) {
        mean_x += x[i];
        mean_y += y[i];
    }

    mean_x /= n;
    mean_y /= n;

    double sum_xx = 0.0;
    double sum_yy = 0.0;
    double sum_xy = 0.0;

    for (std::size_t i = 0; i < x.size(); ++i) {
        const double dx = x[i] - mean_x;
        const double dy = y[i] - mean_y;

        sum_xx += dx * dx;
        sum_yy += dy * dy;
        sum_xy += dx * dy;
    }

    if (sum_xx <= 0.0 || sum_yy <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return sum_xy / std::sqrt(sum_xx * sum_yy);
}


/**
 * Linear interpolation of a sorted time series at time `at`.
 * Clamps outside the observed range. Handles duplicate times.
 */
bool interpolate(
    const double* time,
    const double* value,
    std::size_t n,
    double at,
    double& out
)
{
    if (n == 0) {
        return false;
    }

    if (at <= time[0]) {
        out = value[0];
        return true;
    }

    if (at >= time[n - 1]) {
        out = value[n - 1];
        return true;
    }

    std::size_t lo = 0;
    std::size_t hi = n - 1;

    while (hi - lo > 1) {
        const std::size_t mid = (lo + hi) / 2;

        if (time[mid] <= at) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    if (time[hi] == time[lo]) {
        out = value[lo];
        return true;
    }

    const double frac =
        (at - time[lo]) / (time[hi] - time[lo]);

    out = value[lo] + frac * (value[hi] - value[lo]);

    return true;
}


/**
 * ICCF: Pearson correlation of the continuum with the response
 * interpolated at (continuum_time + lag).
 */
double iccf_value(
    const agnsf::LightCurve& continuum,
    const agnsf::LightCurve& response,
    double lag,
    std::size_t min_overlap,
    std::size_t& overlap
)
{
    const std::size_t nc = continuum.size();
    const std::size_t nr = response.size();

    std::vector<double> x;
    std::vector<double> y;

    x.reserve(nc);
    y.reserve(nc);

    for (std::size_t i = 0; i < nc; ++i) {

        double response_value = 0.0;

        if (!interpolate(
                response.time_data(),
                response.value_data(),
                nr,
                continuum.time_data()[i] + lag,
                response_value
            )) {
            continue;
        }

        x.push_back(continuum.value_data()[i]);
        y.push_back(response_value);
    }

    overlap = x.size();

    if (overlap < min_overlap) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    return pearson(x, y);
}


/**
 * DCF: bin every continuum-response pair by its lag and average the
 * unbinned correlation statistic in each bin.
 */
void dcf_values(
    const agnsf::LightCurve& continuum,
    const agnsf::LightCurve& response,
    const LagGrid& grid,
    double bin_width,
    std::vector<double>& ccf,
    std::vector<std::size_t>& count
)
{
    const std::size_t nc = continuum.size();
    const std::size_t nr = response.size();
    const std::size_t n_taus = ccf.size();

    // Means and standard deviations of the two curves.
    double mean_c = 0.0;
    double mean_r = 0.0;

    for (std::size_t i = 0; i < nc; ++i) mean_c += continuum.value_data()[i];
    for (std::size_t i = 0; i < nr; ++i) mean_r += response.value_data()[i];

    mean_c /= static_cast<double>(nc);
    mean_r /= static_cast<double>(nr);

    double var_c = 0.0;
    double var_r = 0.0;

    for (std::size_t i = 0; i < nc; ++i) {
        const double d = continuum.value_data()[i] - mean_c;
        var_c += d * d;
    }

    for (std::size_t i = 0; i < nr; ++i) {
        const double d = response.value_data()[i] - mean_r;
        var_r += d * d;
    }

    const double sd_c = std::sqrt(var_c / static_cast<double>(nc));
    const double sd_r = std::sqrt(var_r / static_cast<double>(nr));

    if (sd_c == 0.0 || sd_r == 0.0) {
        return; // degenerate curves -> all bins stay NaN
    }

    std::vector<double> sum_udcf(n_taus, 0.0);
    std::vector<double> sum_udcf2(n_taus, 0.0);
    count.assign(n_taus, 0);

    for (std::size_t i = 0; i < nc; ++i) {
        for (std::size_t j = 0; j < nr; ++j) {

            const double lag =
                response.time_data()[j] - continuum.time_data()[i];

            // Bin centered on each grid lag.
            const double k_bin =
                std::round((lag - grid.min) / grid.step);

            if (k_bin < 0.0 ||
                k_bin >= static_cast<double>(n_taus)) {
                continue;
            }

            // |lag - tau_k| must be within half a bin width.
            const std::size_t k = static_cast<std::size_t>(k_bin);

            const double tau_k =
                grid.min + static_cast<double>(k) * grid.step;

            if (std::abs(lag - tau_k) > bin_width / 2.0) {
                continue;
            }

            const double udcf =
                (
                    (continuum.value_data()[i] - mean_c) *
                    (response.value_data()[j] - mean_r)
                ) / (sd_c * sd_r);

            sum_udcf[k] += udcf;
            sum_udcf2[k] += udcf * udcf;
            ++count[k];
        }
    }

    for (std::size_t k = 0; k < n_taus; ++k) {
        if (count[k] == 0) {
            ccf[k] = std::numeric_limits<double>::quiet_NaN();
        } else {
            ccf[k] = sum_udcf[k] / static_cast<double>(count[k]);
        }
    }
}


void estimate_lags(
    LagResult& result,
    double centroid_threshold
)
{
    // Peak.
    double best = -std::numeric_limits<double>::infinity();

    for (std::size_t k = 0; k < result.tau.size(); ++k) {
        if (is_finite(result.ccf[k]) && result.ccf[k] > best) {
            best = result.ccf[k];
            result.lag_peak = result.tau[k];
            result.peak_value = best;
        }
    }

    if (!is_finite(result.peak_value)) {
        return;
    }

    // Centroid above the threshold.
    double sum_w = 0.0;
    double sum_wt = 0.0;

    for (std::size_t k = 0; k < result.tau.size(); ++k) {
        if (is_finite(result.ccf[k]) &&
            result.ccf[k] >= centroid_threshold * result.peak_value) {

            sum_w += result.ccf[k];
            sum_wt += result.tau[k] * result.ccf[k];
        }
    }

    if (sum_w > 0.0) {
        result.lag_centroid = sum_wt / sum_w;
    }
}

} // namespace


std::size_t LagGrid::size() const noexcept
{
    return LagGrid_size(*this);
}


LagResult cross_correlate(
    const agnsf::LightCurve& continuum,
    const agnsf::LightCurve& response,
    const LagGrid& grid,
    const CrossCorrelationConfig& config
)
{
    if (LagGrid_size(grid) == 0) {
        throw std::invalid_argument(
            "timedelay: invalid lag grid (need step > 0, max > min)"
        );
    }

    if (continuum.size() == 0 || response.size() == 0) {
        throw std::invalid_argument(
            "timedelay: light curves must not be empty"
        );
    }

    if (config.dcf_bin_width <= 0.0 ||
        config.centroid_threshold < 0.0 ||
        config.centroid_threshold > 1.0) {
        throw std::invalid_argument(
            "timedelay: invalid cross-correlation config"
        );
    }

    LagResult result;
    result.tau = make_taus(grid);
    result.ccf.assign(result.tau.size(), 0.0);
    result.count.assign(result.tau.size(), 0);

    switch (config.method) {

        case CrossCorrelationMethod::Dcf:
            dcf_values(
                continuum,
                response,
                grid,
                config.dcf_bin_width,
                result.ccf,
                result.count
            );
            break;

        case CrossCorrelationMethod::Iccf:
            for (std::size_t k = 0; k < result.tau.size(); ++k) {
                result.ccf[k] = iccf_value(
                    continuum,
                    response,
                    result.tau[k],
                    config.min_overlap,
                    result.count[k]
                );
            }
            break;
    }

    estimate_lags(result, config.centroid_threshold);

    return result;
}


double lag_value(
    const LagResult& result,
    LagEstimate estimate
)
{
    return estimate == LagEstimate::Centroid
        ? result.lag_centroid
        : result.lag_peak;
}

} // namespace timedelay
} // namespace agnsf
