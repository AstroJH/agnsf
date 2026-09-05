#include <variability/variability.hpp>

#include <cmath>
#include <limits>

namespace agnsf {
namespace variability {

namespace {

double nan()
{
    return std::numeric_limits<double>::quiet_NaN();
}


void set_symmetric(Uncertainty& u, double value, double error)
{
    if (error == error) { // NaN check
        u.lower = value - error;
        u.upper = value + error;
    }
}


/**
 * Regularized upper incomplete gamma function Q(a, x):
 *
 *   Q(a, x) = Gamma(a, x) / Gamma(a)
 *
 * evaluated with the standard series (x < a + 1) and continued
 * fraction (x >= a + 1) expansions (Numerical Recipes, gammp/gammq).
 * Used for the chi^2 variability significance:
 *
 *   Q = Q((n - 1)/2, chi^2/2)
 */
double regularized_upper_incomplete_gamma(double a, double x)
{
    if (!(a > 0.0) || !(x >= 0.0)) {
        return nan();
    }

    if (x == 0.0) {
        return 1.0;
    }

    constexpr int kMaxIterations = 200;
    constexpr double kEps = 3e-14;
    constexpr double kFpMin = 1e-300;

    // Series expansion of P(a, x) (valid for x < a + 1).
    auto series = [&]() {
        double ap = a;
        double sum = 1.0 / a;
        double del = sum;

        for (int i = 1; i <= kMaxIterations; ++i) {
            ap += 1.0;
            del *= x / ap;
            sum += del;

            if (std::fabs(del) < std::fabs(sum) * kEps) {
                break;
            }
        }

        return sum *
            std::exp(-x + a * std::log(x) - std::lgamma(a));
    };

    // Continued fraction for Q(a, x) (valid for x >= a + 1).
    auto continued_fraction = [&]() {
        double b = x + 1.0 - a;
        double c = 1.0 / kFpMin;
        double d = 1.0 / b;
        double h = d;

        for (int i = 1; i <= kMaxIterations; ++i) {
            const double an = -static_cast<double>(i) *
                (static_cast<double>(i) - a);

            b += 2.0;
            d = an * d + b;

            if (std::fabs(d) < kFpMin) {
                d = kFpMin;
            }

            c = b + an / c;

            if (std::fabs(c) < kFpMin) {
                c = kFpMin;
            }

            d = 1.0 / d;

            const double del = d * c;
            h *= del;

            if (std::fabs(del - 1.0) < kEps) {
                break;
            }
        }

        return std::exp(-x + a * std::log(x) - std::lgamma(a)) * h;
    };

    if (x < a + 1.0) {
        return 1.0 - series();
    }

    return continued_fraction();
}

} // namespace


Statistics measure(
    const LightCurveView& curve,
    const Options& options
)
{
    Statistics out;

    const std::size_t n = curve.size();
    out.n = n;

    if (n < 2) {
        return out; // valid == false
    }

    const double* value = curve.value_data();
    const double* error = curve.error_data();

    // Mean, range, mean square error.
    double sum = 0.0;
    double min_value = value[0];
    double max_value = value[0];
    double sum_error2 = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        sum += value[i];
        min_value = std::min(min_value, value[i]);
        max_value = std::max(max_value, value[i]);
        sum_error2 += error[i] * error[i];
    }

    const double err_sys = options.err_sys;
    const double mean = sum / static_cast<double>(n);
    const double mean_square_error =
        sum_error2 / static_cast<double>(n) +
        err_sys * err_sys;

    // Variance (sample, ddof = 1) around the chosen mean.
    double mean_for_variance = mean;

    if (options.weighted) {
        // Inverse-variance weighted mean; points with sigma <= 0 are
        // skipped (no information).
        double sum_w = 0.0;
        double sum_wx = 0.0;

        for (std::size_t i = 0; i < n; ++i) {
            if (error[i] > 0.0) {
                const double w = 1.0 / (error[i] * error[i]);
                sum_w += w;
                sum_wx += w * value[i];
            }
        }

        if (sum_w > 0.0) {
            out.weighted_mean = sum_wx / sum_w;
            out.weighted_mean_error = 1.0 / std::sqrt(sum_w);
            mean_for_variance = out.weighted_mean;
        } else {
            out.weighted_mean = nan();
            out.weighted_mean_error = nan();
        }
    } else {
        out.weighted_mean = nan();
        out.weighted_mean_error = nan();
    }

    double sum_square = 0.0;
    double sum_successive_square = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        const double diff = value[i] - mean_for_variance;
        sum_square += diff * diff;

        if (i + 1 < n) {
            const double step = value[i + 1] - value[i];
            sum_successive_square += step * step;
        }
    }

    const double sample_variance =
        sum_square / static_cast<double>(n - 1);

    const double stddev = std::sqrt(sample_variance);

    out.mean = mean;
    out.stddev = stddev;
    out.stddev_error =
        stddev / std::sqrt(2.0 * static_cast<double>(n - 1));

    out.peak_to_peak = max_value - min_value;

    const double range_square =
        (max_value - min_value) * (max_value - min_value);

    out.peak_to_peak_noise_corrected = range_square > 2.0 * mean_square_error
        ? std::sqrt(range_square - 2.0 * mean_square_error)
        : 0.0;

    // Intrinsic amplitude and excess variance.
    const double excess = sample_variance - mean_square_error;

    out.sigma_m = excess > 0.0 ? std::sqrt(excess) : 0.0;

    out.xs = excess;

    out.nxs = mean != 0.0 ? excess / (mean * mean) : nan();
    out.fvar = (out.nxs == out.nxs && out.nxs > 0.0)
        ? std::sqrt(out.nxs)
        : 0.0;

    // Vaughan et al. (2003) analytic uncertainties.
    if (mean != 0.0 && out.fvar > 0.0) {
        const double nn = static_cast<double>(n);
        const double mean2 = mean * mean;

        const double fvar_error = std::sqrt(
            mean_square_error / (nn * mean2) +
            0.5 / nn *
                std::pow(mean_square_error / (mean2 * out.fvar), 2.0)
        );

        set_symmetric(
            out.fvar_uncertainty, out.fvar, fvar_error
        );

        const double nxs_error = 2.0 * out.fvar * fvar_error;

        set_symmetric(
            out.nxs_uncertainty, out.nxs, nxs_error
        );

        set_symmetric(
            out.xs_uncertainty, out.xs, mean2 * nxs_error
        );
    }

    // chi^2 variability significance.
    std::size_t n_valid_chi2 = 0;
    double chi2 = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        if (error[i] > 0.0) {
            const double variance =
                error[i] * error[i] +
                options.err_sys * options.err_sys;

            const double residual = value[i] - mean;
            chi2 += residual * residual / variance;
            ++n_valid_chi2;
        }
    }

    if (n_valid_chi2 >= 2) {
        const double dof = static_cast<double>(n_valid_chi2 - 1);

        out.chi2 = chi2;
        out.chi2_dof = dof;
        out.chi2_q = regularized_upper_incomplete_gamma(
            dof / 2.0,
            chi2 / 2.0
        );
    } else {
        out.chi2 = nan();
        out.chi2_dof = nan();
        out.chi2_q = nan();
    }

    // Von Neumann ratio: mean square successive difference / variance
    // (only meaningful for time-ordered, evenly sampled data).
    if (sum_square > 0.0) {
        out.von_neumann = sum_successive_square / sum_square;
    } else {
        out.von_neumann = nan();
    }

    out.valid = true;
    return out;
}

} // namespace variability
} // namespace agnsf
