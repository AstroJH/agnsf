#include <timedelay/transfer_function.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <vector>

namespace agnsf {
namespace timedelay {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kGaussianSupport = 4.0; // sigma multiples

double gaussian_support_width(double width)
{
    return kGaussianSupport * width;
}


double psi(
    double tau,
    TransferFunctionShape shape,
    double lag,
    double width
)
{
    const double z = tau - lag;

    if (shape == TransferFunctionShape::TopHat) {
        if (std::abs(z) > width) {
            return 0.0;
        }

        return 1.0 / (2.0 * width);
    }

    // Gaussian
    const double t = z / width;

    return std::exp(-0.5 * t * t) /
        (std::sqrt(2.0 * kPi) * width);
}


void validate_finite(const agnsf::LightCurve& curve)
{
    const std::size_t n = curve.size();

    for (std::size_t i = 0; i < n; ++i) {
        if (!std::isfinite(curve.time_data()[i]) ||
            !std::isfinite(curve.value_data()[i]) ||
            !std::isfinite(curve.error_data()[i])) {

            throw std::invalid_argument(
                "transfer function: light curves must "
                "contain finite data"
            );
        }
    }
}


/**
 * Uniform grid over the continuum time range with linearly
 * interpolated continuum values. The grid is built once and reused
 * by every objective evaluation, so the optimization loop performs
 * no per-evaluation interpolation setup.
 */
struct ContinuumGrid {
    std::vector<double> s;
    std::vector<double> c;
    double step = 0.0;
};


double median_spacing(
    const agnsf::LightCurve& continuum
)
{
    const std::size_t n = continuum.size();

    std::vector<double> spacings;
    spacings.reserve(n - 1);

    for (std::size_t i = 1; i < n; ++i) {
        spacings.push_back(
            continuum.time_data()[i] - continuum.time_data()[i - 1]
        );
    }

    std::sort(spacings.begin(), spacings.end());

    const std::size_t mid = spacings.size() / 2;

    if (spacings.size() % 2 == 0) {
        return 0.5 * (spacings[mid - 1] + spacings[mid]);
    }

    return spacings[mid];
}


ContinuumGrid build_continuum_grid(
    const agnsf::LightCurve& continuum,
    double grid_step
)
{
    const std::size_t n = continuum.size();

    const double t_min = continuum.time_data()[0];
    const double t_max = continuum.time_data()[n - 1];

    if (t_max <= t_min) {
        throw std::invalid_argument(
            "transfer function: continuum needs a "
            "non-degenerate time range"
        );
    }

    double step = grid_step;

    if (step <= 0.0) {
        // Default: a few samples per median continuum cadence.
        step = median_spacing(continuum) / 5.0;
    }

    if (step <= 0.0 || !std::isfinite(step)) {
        throw std::invalid_argument(
            "transfer function: invalid continuum grid step"
        );
    }

    const std::size_t count = static_cast<std::size_t>(
        std::floor((t_max - t_min) / step)
    ) + 1;

    ContinuumGrid grid;
    grid.step = step;
    grid.s.reserve(count);
    grid.c.reserve(count);

    const double* time = continuum.time_data();
    const double* value = continuum.value_data();

    for (std::size_t i = 0; i < count; ++i) {
        const double s = t_min + static_cast<double>(i) * step;

        grid.s.push_back(s);

        // Linear interpolation between bracketing samples.
        const auto it = std::lower_bound(
            time,
            time + n,
            s
        );

        if (it == time) {
            grid.c.push_back(value[0]);
        } else if (it == time + n) {
            grid.c.push_back(value[n - 1]);
        } else {
            const std::size_t hi =
                static_cast<std::size_t>(it - time);
            const std::size_t lo = hi - 1;

            const double t_lo = time[lo];
            const double t_hi = time[hi];
            const double f = (s - t_lo) / (t_hi - t_lo);

            grid.c.push_back(
                value[lo] + f * (value[hi] - value[lo])
            );
        }
    }

    return grid;
}


/**
 * Model response at one time:
 *
 *   R(t) = offset + amplitude * int [ Psi(tau) C(t - tau) ] d tau
 *
 * The integral is approximated on the fixed continuum grid with the
 * rectangle rule over the transfer-function support
 * [lag - k*width, lag + k*width]. The continuum is extended flat
 * outside its observed range, so every response point always has a
 * well-defined model value and ALL response points contribute to
 * chi^2. (At the optimum the windows lie inside the observed range,
 * where the extension is never used.)
 */
double model_response_at(
    double t,
    const ContinuumGrid& grid,
    TransferFunctionShape shape,
    double offset,
    double amplitude,
    double lag,
    double width
)
{
    const double half =
        shape == TransferFunctionShape::TopHat
            ? width
            : gaussian_support_width(width);

    const double s_lo = t - (lag + half);
    const double s_hi = t - (lag - half);

    const double s_min = grid.s.front();
    const double s_max = grid.s.back();

    // Integer grid range covering [s_lo, s_hi]; indices outside
    // [0, size) use the flat extension of the continuum.
    const long long i_lo = static_cast<long long>(
        std::floor((s_lo - s_min) / grid.step)
    );
    const long long i_hi = static_cast<long long>(
        std::ceil((s_hi - s_min) / grid.step)
    );

    double integral = 0.0;

    for (long long i = i_lo; i <= i_hi; ++i) {
        double c;

        if (i < 0) {
            c = grid.c.front();
        } else if (static_cast<std::size_t>(i) >= grid.c.size()) {
            c = grid.c.back();
        } else {
            c = grid.c[static_cast<std::size_t>(i)];
        }

        const double s = s_min + static_cast<double>(i) * grid.step;
        const double tau = t - s;

        integral += c * psi(tau, shape, lag, width);
    }

    integral *= grid.step;

    return offset + amplitude * integral;
}


std::size_t count_covered_points(
    const agnsf::LightCurve& response,
    const ContinuumGrid& grid,
    TransferFunctionShape shape,
    double lag,
    double width
)
{
    const double half =
        shape == TransferFunctionShape::TopHat
            ? width
            : gaussian_support_width(width);

    const double s_min = grid.s.front();
    const double s_max = grid.s.back();

    std::size_t count = 0;

    for (std::size_t j = 0; j < response.size(); ++j) {
        const double t = response.time_data()[j];

        const double window_lo = t - (lag + half);
        const double window_hi = t - (lag - half);

        if (window_lo >= s_min && window_hi <= s_max) {
            ++count;
        }
    }

    return count;
}


/**
 * Transfer-function objective with the linear parameters profiled.
 *
 * The model is linear in (offset, amplitude) for fixed (lag, width):
 *
 *   R_j = offset + amplitude * k_j,   k_j = (C * Psi)(t_j)
 *
 * so the best (offset, amplitude) for a given (lag, width) is the
 * weighted linear least-squares solution. Profiling removes the two
 * linear parameters from the non-linear search, which makes the
 * (lag, width) surface much more regular and the fit far less
 * sensitive to starting values. The solver (BOBYQA) only varies the
 * free non-linear parameters.
 *
 * Assumption (kept for simplicity, documented in the header):
 * measurement contributions of different pairs are independent and
 * only the response errors enter chi^2 (the continuum is treated as
 * known).
 */
class ProfiledTransferFunctionObjective {
public:
    ProfiledTransferFunctionObjective(
        ContinuumGrid grid,
        const agnsf::LightCurve& response,
        TransferFunctionShape shape,
        bool offset_free,
        double offset_pinned,
        double offset_lower,
        double offset_upper,
        bool amplitude_free,
        double amplitude_pinned,
        double amplitude_lower,
        double amplitude_upper
    )
        : grid_(std::move(grid)),
          shape_(shape),
          offset_free_(offset_free),
          offset_pinned_(offset_pinned),
          offset_lower_(offset_lower),
          offset_upper_(offset_upper),
          amplitude_free_(amplitude_free),
          amplitude_pinned_(amplitude_pinned),
          amplitude_lower_(amplitude_lower),
          amplitude_upper_(amplitude_upper)
    {
        times_.reserve(response.size());
        values_.reserve(response.size());
        errors_.reserve(response.size());

        for (std::size_t i = 0; i < response.size(); ++i) {
            times_.push_back(response.time_data()[i]);
            values_.push_back(response.value_data()[i]);
            errors_.push_back(response.error_data()[i]);
        }

        k_.resize(response.size());
    }

    /**
     * chi^2 at (lag, width) with (offset, amplitude) profiled (or
     * pinned); optionally returns the profiled linear parameters.
     */
    double solve(
        double lag,
        double width,
        double* offset_out,
        double* amplitude_out
    ) const
    {
        // First pass: convolution kernels k_j and weighted sums.
        double sum_w = 0.0;
        double sum_wk = 0.0;
        double sum_wkk = 0.0;
        double sum_wr = 0.0;
        double sum_wkr = 0.0;

        for (std::size_t j = 0; j < times_.size(); ++j) {
            const double k = model_response_at(
                times_[j],
                grid_,
                shape_,
                0.0,
                1.0,
                lag,
                width
            );

            k_[j] = k;

            const double sigma = errors_[j];
            const double weight =
                sigma > 0.0 ? 1.0 / (sigma * sigma) : 1.0;

            const double r = values_[j];

            sum_w += weight;
            sum_wk += weight * k;
            sum_wkk += weight * k * k;
            sum_wr += weight * r;
            sum_wkr += weight * k * r;
        }

        double offset = offset_pinned_;
        double amplitude = amplitude_pinned_;

        if (offset_free_ && amplitude_free_) {
            const double det = sum_w * sum_wkk - sum_wk * sum_wk;

            if (std::abs(det) >
                1e-12 * std::max(sum_w * sum_wkk, 1.0)) {

                offset = (sum_wkk * sum_wr - sum_wk * sum_wkr) / det;
                amplitude = (sum_w * sum_wkr - sum_wk * sum_wr) / det;
            } else {
                // k is (nearly) constant: amplitude is degenerate.
                offset = sum_w > 0.0 ? sum_wr / sum_w : 0.0;
                amplitude = 0.0;
            }
        } else if (offset_free_) {
            offset = sum_w > 0.0
                ? (sum_wr - amplitude * sum_wk) / sum_w
                : 0.0;
        } else if (amplitude_free_) {
            amplitude = sum_wkk > 0.0
                ? (sum_wkr - offset * sum_wk) / sum_wkk
                : 0.0;
        }

        // Respect the user's box bounds on the linear parameters.
        offset = std::clamp(offset, offset_lower_, offset_upper_);
        amplitude = std::clamp(
            amplitude, amplitude_lower_, amplitude_upper_
        );

        // Second pass: chi^2 with the chosen linear parameters.
        double chi2 = 0.0;

        for (std::size_t j = 0; j < times_.size(); ++j) {
            const double sigma = errors_[j];
            const double weight =
                sigma > 0.0 ? 1.0 / (sigma * sigma) : 1.0;

            const double residual =
                values_[j] - (offset + amplitude * k_[j]);

            const double z = residual * std::sqrt(weight);
            chi2 += z * z;
        }

        if (offset_out != nullptr) {
            *offset_out = offset;
        }

        if (amplitude_out != nullptr) {
            *amplitude_out = amplitude;
        }

        return chi2;
    }

private:
    ContinuumGrid grid_;
    TransferFunctionShape shape_;

    bool offset_free_;
    double offset_pinned_;
    double offset_lower_;
    double offset_upper_;

    bool amplitude_free_;
    double amplitude_pinned_;
    double amplitude_lower_;
    double amplitude_upper_;

    std::vector<double> times_;
    std::vector<double> values_;
    std::vector<double> errors_;

    mutable std::vector<double> k_;
};


/**
 * Adapts the profiled objective to the generic optimizer over the
 * free non-linear parameters [lag, width] (a subset may be pinned).
 */
class ReducedObjective : public optimization::Objective {
public:
    ReducedObjective(
        const ProfiledTransferFunctionObjective& inner,
        bool lag_free,
        bool width_free,
        double lag_pinned,
        double width_pinned
    )
        : inner_(inner),
          lag_free_(lag_free),
          width_free_(width_free),
          lag_pinned_(lag_pinned),
          width_pinned_(width_pinned)
    {
    }

    double evaluate(
        const std::vector<double>& x
    ) const override
    {
        std::size_t index = 0;

        const double lag =
            lag_free_ ? x[index++] : lag_pinned_;
        const double width =
            width_free_ ? x[index] : width_pinned_;

        return inner_.solve(lag, width, nullptr, nullptr);
    }

private:
    const ProfiledTransferFunctionObjective& inner_;
    bool lag_free_;
    bool width_free_;
    double lag_pinned_;
    double width_pinned_;
};


bool is_pinned(const optimization::Parameter& p)
{
    return p.fixed || p.lower == p.upper;
}

} // namespace


std::vector<optimization::Parameter>
default_transfer_function_parameters(
    const LightCurve& continuum,
    const LightCurve& response,
    const TransferFunctionConfig& config
)
{
    if (continuum.size() < 2 || response.size() < 1) {
        throw std::invalid_argument(
            "transfer function: need at least two continuum "
            "samples and one response sample"
        );
    }

    validate_finite(continuum);
    validate_finite(response);

    const double* tc = continuum.time_data();
    const double* vc = continuum.value_data();
    const double* tr = response.time_data();
    const double* vr = response.value_data();

    const std::size_t nc = continuum.size();
    const std::size_t nr = response.size();

    // Physically meaningful lag range from the data windows.
    double lag_lo = tr[0] - tc[nc - 1];
    double lag_hi = tr[nr - 1] - tc[0];

    if (lag_hi < lag_lo) {
        std::swap(lag_lo, lag_hi);
    }

    if (lag_hi == lag_lo) {
        lag_lo -= 1.0;
        lag_hi += 1.0;
    }

    double min_c = vc[0];
    double max_c = vc[0];
    double min_r = vr[0];
    double max_r = vr[0];
    double sum_c = 0.0;
    double sum_r = 0.0;

    for (std::size_t i = 0; i < nc; ++i) {
        min_c = std::min(min_c, vc[i]);
        max_c = std::max(max_c, vc[i]);
        sum_c += vc[i];
    }

    for (std::size_t i = 0; i < nr; ++i) {
        min_r = std::min(min_r, vr[i]);
        max_r = std::max(max_r, vr[i]);
        sum_r += vr[i];
    }

    const double range_c = std::max(max_c - min_c, 1e-12);
    const double range_r = std::max(max_r - min_r, 1e-12);

    const double amplitude0 = range_r / range_c;
    const double mean_c = sum_c / static_cast<double>(nc);
    const double mean_r = sum_r / static_cast<double>(nr);

    // Grid step (needed for a positive width bound).
    ContinuumGrid grid =
        build_continuum_grid(continuum, config.grid_step);

    // Width defaults are tied to the continuum cadence: a transfer
    // function much broader than the typical variability timescale is
    // usually unphysical, and starting from a very broad width traps
    // the local optimizer far from the best fit.
    const double median = median_spacing(continuum);

    const double width_lo = grid.step;
    const double width_hi = std::max(
        10.0 * median,
        0.5 * (lag_hi - lag_lo)
    );
    const double width_init = std::clamp(
        5.0 * median,
        width_lo,
        width_hi
    );

    std::vector<optimization::Parameter> parameters(4);

    // [0] offset
    parameters[0].value = mean_r - amplitude0 * mean_c;
    parameters[0].lower = min_r - 2.0 * range_r;
    parameters[0].upper = max_r + 2.0 * range_r;

    // [1] amplitude
    parameters[1].value = amplitude0;
    parameters[1].lower = 0.0;
    parameters[1].upper = 100.0 * amplitude0;

    // [2] lag
    parameters[2].value = 0.5 * (lag_lo + lag_hi);
    parameters[2].lower = lag_lo;
    parameters[2].upper = lag_hi;

    // [3] width
    parameters[3].value = width_init;
    parameters[3].lower = width_lo;
    parameters[3].upper = width_hi;

    return parameters;
}


TransferFunctionResult fit_transfer_function(
    const LightCurve& continuum,
    const LightCurve& response,
    const std::vector<optimization::Parameter>& parameters,
    const TransferFunctionConfig& config,
    const optimization::Options& options
)
{
    if (continuum.size() < 2 || response.size() < 1) {
        throw std::invalid_argument(
            "transfer function: need at least two continuum "
            "samples and one response sample"
        );
    }

    validate_finite(continuum);
    validate_finite(response);

    if (parameters.size() != kTransferFunctionParameters) {
        throw std::invalid_argument(
            "transfer function: expected 4 parameters "
            "[offset, amplitude, lag, width]"
        );
    }

    if (config.grid_step < 0.0 ||
        !std::isfinite(config.grid_step)) {

        throw std::invalid_argument(
            "transfer function: grid_step must be >= 0"
        );
    }

    for (const auto& p : parameters) {
        if (!std::isfinite(p.value) ||
            !std::isfinite(p.lower) ||
            !std::isfinite(p.upper)) {

            throw std::invalid_argument(
                "transfer function: parameter values and "
                "bounds must be finite"
            );
        }

        if (p.lower > p.upper) {
            throw std::invalid_argument(
                "transfer function: lower bound exceeds "
                "upper bound"
            );
        }
    }

    if (parameters[3].lower <= 0.0) {
        throw std::invalid_argument(
            "transfer function: width lower bound must "
            "be positive"
        );
    }

    const ContinuumGrid grid =
        build_continuum_grid(continuum, config.grid_step);

    const bool offset_free = !is_pinned(parameters[0]);
    const bool amplitude_free = !is_pinned(parameters[1]);
    const bool lag_free = !is_pinned(parameters[2]);
    const bool width_free = !is_pinned(parameters[3]);

    const ProfiledTransferFunctionObjective inner(
        grid,
        response,
        config.shape,
        offset_free,
        parameters[0].value,
        parameters[0].lower,
        parameters[0].upper,
        amplitude_free,
        parameters[1].value,
        parameters[1].lower,
        parameters[1].upper
    );

    // Free non-linear parameters passed to the generic optimizer.
    std::vector<optimization::Parameter> nonlinear;

    if (lag_free) {
        nonlinear.push_back(parameters[2]);
    }

    if (width_free) {
        nonlinear.push_back(parameters[3]);
    }

    const ReducedObjective objective(
        inner,
        lag_free,
        width_free,
        parameters[2].value,
        parameters[3].value
    );

    if (nonlinear.empty()) {
        // All non-linear parameters are pinned: a single evaluation
        // with the (profiled or pinned) linear parameters.
        TransferFunctionResult result;

        result.converged = true;
        result.message = "non-linear parameters fixed";
        result.evaluations = 1;
        result.lag = parameters[2].value;
        result.width = parameters[3].value;

        result.chi2 = inner.solve(
            result.lag,
            result.width,
            &result.offset,
            &result.amplitude
        );

        result.n_valid_points = count_covered_points(
            response,
            grid,
            config.shape,
            result.lag,
            result.width
        );

        return result;
    }

    optimization::Result best;
    bool have_best = false;

    const std::size_t restarts = config.lag_restarts;

    if (lag_free && restarts != 1) {
        // Coarse-to-fine lag scan: rank lag grid points by the
        // profiled chi^2 (with the initial width), then refine the
        // best candidates with the local optimizer. The chi^2 vs lag
        // surface is multimodal for quasi-periodic continua, so a
        // single local optimization is not reliable.
        const double lag_lo = parameters[2].lower;
        const double lag_hi = parameters[2].upper;

        const double width_init = std::clamp(
            parameters[3].value,
            parameters[3].lower,
            parameters[3].upper
        );

        std::size_t count = restarts;

        if (count == 0) {
            const double spacing =
                std::max(2.0 * width_init, grid.step);

            count = static_cast<std::size_t>(
                std::ceil((lag_hi - lag_lo) / spacing)
            );
            count = std::clamp<std::size_t>(count, 8, 64);
        }

        if (count < 2) {
            count = 2;
        }

        std::vector<std::pair<double, double>> scored; // (chi2, lag)

        for (std::size_t i = 0; i < count; ++i) {
            const double fraction =
                static_cast<double>(i) /
                static_cast<double>(count - 1);

            const double candidate =
                lag_lo + fraction * (lag_hi - lag_lo);

            const double chi2 =
                inner.solve(candidate, width_init, nullptr, nullptr);

            scored.emplace_back(chi2, candidate);
        }

        std::sort(scored.begin(), scored.end());

        const std::size_t refine =
            std::min<std::size_t>(scored.size(), 3);

        for (std::size_t i = 0; i < refine; ++i) {
            nonlinear[0].value = scored[i].second;

            const optimization::Result result =
                optimization::minimize(objective, nonlinear, options);

            if (!have_best ||
                result.objective_value < best.objective_value) {

                best = result;
                have_best = true;
            }
        }
    } else {
        const optimization::Result result =
            optimization::minimize(objective, nonlinear, options);

        best = result;
        have_best = true;
    }

    TransferFunctionResult result;

    result.converged = have_best && best.converged;
    result.message = best.message;
    result.evaluations = best.evaluations;

    if (have_best) {
        std::size_t index = 0;

        const double lag =
            lag_free ? best.parameters[index++] : parameters[2].value;
        const double width =
            width_free ? best.parameters[index] : parameters[3].value;

        result.lag = lag;
        result.width = width;

        result.chi2 = inner.solve(
            lag,
            width,
            &result.offset,
            &result.amplitude
        );

        result.n_valid_points = count_covered_points(
            response,
            grid,
            config.shape,
            lag,
            width
        );
    }

    return result;
}


double evaluate_transfer_function(
    double tau,
    TransferFunctionShape shape,
    double lag,
    double width
)
{
    if (!std::isfinite(tau) ||
        !std::isfinite(lag) ||
        !std::isfinite(width) ||
        width <= 0.0) {

        throw std::invalid_argument(
            "transfer function: tau, lag and width must be "
            "finite with width > 0"
        );
    }

    return psi(tau, shape, lag, width);
}


std::vector<double> transfer_function_curve(
    const std::vector<double>& taus,
    TransferFunctionShape shape,
    double lag,
    double width
)
{
    std::vector<double> curve;
    curve.reserve(taus.size());

    for (const double tau : taus) {
        curve.push_back(
            evaluate_transfer_function(tau, shape, lag, width)
        );
    }

    return curve;
}


std::vector<double> transfer_function_model_response(
    const LightCurve& continuum,
    const std::vector<double>& response_times,
    TransferFunctionShape shape,
    double offset,
    double amplitude,
    double lag,
    double width,
    double grid_step
)
{
    if (continuum.size() < 2) {
        throw std::invalid_argument(
            "transfer function: need at least two continuum "
            "samples"
        );
    }

    validate_finite(continuum);

    if (width <= 0.0 || !std::isfinite(width)) {
        throw std::invalid_argument(
            "transfer function: width must be finite and > 0"
        );
    }

    const ContinuumGrid grid =
        build_continuum_grid(continuum, grid_step);

    std::vector<double> model;
    model.reserve(response_times.size());

    for (const double t : response_times) {
        model.push_back(
            model_response_at(
                t,
                grid,
                shape,
                offset,
                amplitude,
                lag,
                width
            )
        );
    }

    return model;
}

} // namespace timedelay
} // namespace agnsf
