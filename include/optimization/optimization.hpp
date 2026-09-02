#pragma once

#include <cstddef>
#include <string>
#include <vector>

namespace agnsf {
namespace optimization {

/**
 * Bound-constrained minimizer selection.
 *
 * Both algorithms are provided by NLopt; this module only wires a
 * generic Objective to the library.
 */
enum class Algorithm {
    // NLopt LN_BOBYQA: derivative-free, bound-constrained local
    // optimizer based on Powell's quadratic approximation. Default.
    BoundedBobyqa = 0,

    // NLopt LN_NELDERMEAD: derivative-free simplex with bounds.
    NelderMead
};


/**
 * One model parameter: initial value, box bounds, and optional pinning.
 *
 * A parameter is kept fixed when `fixed == true` or when
 * `lower == upper`.
 */
struct Parameter {
    double value = 0.0;
    double lower = 0.0;
    double upper = 0.0;
    bool fixed = false;
};


/**
 * Solver options.
 */
struct Options {
    Algorithm algorithm = Algorithm::BoundedBobyqa;

    // Relative tolerance on the parameters.
    double xtol_rel = 1e-6;

    // Relative tolerance on the objective value.
    double ftol_rel = 1e-8;

    // Maximum number of objective evaluations.
    std::size_t max_evaluations = 10000;
};


enum class Status {
    Success = 0,    // converged within the requested tolerances
    MaxEvaluations, // stopped by the evaluation budget
    Invalid,        // invalid input (bad bounds, size, ...)
    Failure         // the optimizer reported an error
};


/**
 * Result of a minimization.
 *
 * `parameters` always holds the full parameter vector (fixed
 * parameters pinned at their initial value).
 */
struct Result {
    Status status = Status::Failure;
    bool converged = false;
    std::vector<double> parameters;
    double objective_value = 0.0;
    std::size_t evaluations = 0;
    std::string message;
};


/**
 * Scalar objective to minimize over the full parameter vector.
 *
 * `evaluate` receives the full vector (fixed parameters substituted),
 * so domain code can index parameters by position/name without
 * knowing which parameters were pinned.
 *
 * The optimization loop calls only C++ code; no Python callbacks are
 * involved, so a fit never round-trips through the interpreter.
 */
class Objective {
public:
    virtual ~Objective() = default;
    virtual double evaluate(const std::vector<double>& x) const = 0;
};


/**
 * Minimize `objective` over the box-constrained `parameters`.
 *
 * Fixed parameters (fixed == true or lower == upper) are pinned.
 * Initial values are clamped into [lower, upper] before optimizing.
 *
 * @throws std::invalid_argument for invalid input:
 *   - parameters is empty
 *   - lower > upper for any parameter
 *   - any initial value / bound is non-finite
 */
Result minimize(
    const Objective& objective,
    const std::vector<Parameter>& parameters,
    const Options& options = {}
);

} // namespace optimization
} // namespace agnsf
