#pragma once

#include <cstddef>
#include <string>
#include <vector>

#include <core/light_curve.hpp>
#include <optimization/optimization.hpp>

namespace agnsf {
namespace timedelay {

/**
 * Parametric transfer-function shape.
 *
 * Both shapes integrate to unity, so the fitted `amplitude` is the
 * integrated response (gain):
 *
 *   Gaussian:  Psi(tau) = exp(-(tau - lag)^2 / (2 width^2))
 *                          / (sqrt(2 pi) width)
 *
 *   TopHat:    Psi(tau) = 1 / (2 width)  for |tau - lag| <= width
 *                         0              otherwise
 *
 * The model response is
 *
 *   R(t) = offset + amplitude * (C * Psi)(t)
 *
 * where C is the (interpolated) continuum light curve.
 */
enum class TransferFunctionShape {
    Gaussian = 0,
    TopHat
};


struct TransferFunctionConfig {
    TransferFunctionShape shape = TransferFunctionShape::Gaussian;

    // Continuum interpolation step used for the convolution integral.
    // 0 selects an automatic step (median continuum spacing / 5).
    double grid_step = 0.0;

    // Lag scan for the initial search (chi^2 vs lag is multimodal
    // for quasi-periodic continua):
    //   0 = auto: scan the lag bound with a spacing of ~2*width,
    //       refining the best candidates with the local optimizer;
    //   1 = single start at the provided initial lag;
    //   n > 1 = explicit scan size.
    std::size_t lag_restarts = 0;
};


/**
 * Parameter order expected by fit_transfer_function() and returned by
 * default_transfer_function_parameters():
 *
 *   [0] offset, [1] amplitude, [2] lag, [3] width
 */
constexpr std::size_t kTransferFunctionParameters = 4;


struct TransferFunctionResult {
    bool converged = false;

    double offset = 0.0;
    double amplitude = 0.0;
    double lag = 0.0;
    double width = 0.0;

    double chi2 = 0.0;
    std::size_t evaluations = 0;

    // Number of response points whose convolution window lies fully
    // inside the continuum coverage at the best fit (informational;
    // all response points always contribute to chi^2).
    std::size_t n_valid_points = 0;
    std::string message;
};


/**
 * Data-driven initial values and box bounds for the four transfer
 * function parameters (see kTransferFunctionParameters).
 */
std::vector<optimization::Parameter> default_transfer_function_parameters(
    const LightCurve& continuum,
    const LightCurve& response,
    const TransferFunctionConfig& config = {}
);


/**
 * Fit the parametric transfer function by chi^2 minimization.
 *
 *   chi^2 = sum_j [ (R_j - R_model(t_j)) / sigma_j ]^2
 *
 * over all response points. The continuum is taken as known (only
 * the response errors enter chi^2), linearly interpolated onto a
 * uniform grid, and extended flat outside its observed range so that
 * every response point always contributes. The convolution integral
 * is evaluated on that grid inside the transfer-function support.
 * Points with sigma <= 0 are treated as unweighted (weight 1).
 *
 * The linear parameters (offset, amplitude) are profiled analytically
 * (weighted least squares) for each (lag, width); only the free
 * non-linear parameters are searched by the local optimizer, which
 * makes the fit robust to starting values. A lag scan (see
 * TransferFunctionConfig::lag_restarts) reduces the risk of landing
 * in a secondary minimum of the multimodal chi^2 vs lag surface.
 *
 * @throws std::invalid_argument for invalid input (too few data
 *   points, non-finite data, parameters.size() != 4, lower > upper,
 *   width.lower <= 0, ...).
 */
TransferFunctionResult fit_transfer_function(
    const LightCurve& continuum,
    const LightCurve& response,
    const std::vector<optimization::Parameter>& parameters,
    const TransferFunctionConfig& config = {},
    const optimization::Options& options = {}
);


/**
 * Evaluate the normalized transfer function Psi(tau).
 */
double evaluate_transfer_function(
    double tau,
    TransferFunctionShape shape,
    double lag,
    double width
);


/**
 * Evaluate the normalized transfer function on a grid of lags.
 */
std::vector<double> transfer_function_curve(
    const std::vector<double>& taus,
    TransferFunctionShape shape,
    double lag,
    double width
);


/**
 * Model response at arbitrary times for given transfer-function
 * parameters (used to overlay the fitted model on the data).
 */
std::vector<double> transfer_function_model_response(
    const LightCurve& continuum,
    const std::vector<double>& response_times,
    TransferFunctionShape shape,
    double offset,
    double amplitude,
    double lag,
    double width,
    double grid_step = 0.0
);

} // namespace timedelay
} // namespace agnsf
