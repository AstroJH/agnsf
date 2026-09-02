#include <optimization/optimization.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

#include <nlopt.hpp>

namespace agnsf {
namespace optimization {

namespace {

struct NloptData {
    const Objective* objective;
    std::vector<std::size_t> free_indices;
    std::vector<double> full;
};


// NLopt vfunc trampoline: expand the reduced free-parameter vector to
// the full vector and forward to the domain Objective.
double objective_trampoline(
    const std::vector<double>& x,
    std::vector<double>& grad,
    void* data
)
{
    (void)grad; // derivative-free algorithms only

    NloptData* d = static_cast<NloptData*>(data);

    for (std::size_t i = 0; i < d->free_indices.size(); ++i) {
        d->full[d->free_indices[i]] = x[i];
    }

    return d->objective->evaluate(d->full);
}


void validate_parameters(
    const std::vector<Parameter>& parameters
)
{
    if (parameters.empty()) {
        throw std::invalid_argument(
            "optimization: at least one parameter is required"
        );
    }

    for (const auto& p : parameters) {
        if (!std::isfinite(p.value) ||
            !std::isfinite(p.lower) ||
            !std::isfinite(p.upper)) {

            throw std::invalid_argument(
                "optimization: parameter values and bounds "
                "must be finite"
            );
        }

        if (p.lower > p.upper) {
            throw std::invalid_argument(
                "optimization: lower bound exceeds upper bound"
            );
        }
    }
}


Result evaluate_only(
    const Objective& objective,
    const std::vector<Parameter>& parameters
)
{
    std::vector<double> x;
    x.reserve(parameters.size());

    for (const auto& p : parameters) {
        x.push_back(p.value);
    }

    Result result;
    result.status = Status::Success;
    result.converged = true;
    result.parameters = x;
    result.objective_value = objective.evaluate(x);
    result.evaluations = 1;
    result.message = "all parameters fixed";
    return result;
}

} // namespace


Result minimize(
    const Objective& objective,
    const std::vector<Parameter>& parameters,
    const Options& options
)
{
    validate_parameters(parameters);

    // Free parameters: those not explicitly fixed and with a
    // non-degenerate interval. lower == upper pins the value.
    std::vector<std::size_t> free_indices;

    for (std::size_t i = 0; i < parameters.size(); ++i) {
        const auto& p = parameters[i];

        if (p.fixed || p.lower == p.upper) {
            continue;
        }

        free_indices.push_back(i);
    }

    if (free_indices.empty()) {
        return evaluate_only(objective, parameters);
    }

    NloptData data;
    data.objective = &objective;
    data.free_indices = free_indices;

    // Seed the full vector with the initial values so that fixed
    // parameters keep their pinned value even before the first
    // objective evaluation.
    data.full.resize(parameters.size());

    for (std::size_t i = 0; i < parameters.size(); ++i) {
        data.full[i] = parameters[i].value;
    }

    std::vector<double> lower;
    std::vector<double> upper;
    std::vector<double> x;

    lower.reserve(free_indices.size());
    upper.reserve(free_indices.size());
    x.reserve(free_indices.size());

    for (const std::size_t index : free_indices) {
        const auto& p = parameters[index];

        // Clamp the initial value into the feasible box.
        const double value =
            std::clamp(p.value, p.lower, p.upper);

        lower.push_back(p.lower);
        upper.push_back(p.upper);
        x.push_back(value);
    }

    nlopt::algorithm algorithm =
        options.algorithm == Algorithm::NelderMead
            ? nlopt::LN_NELDERMEAD
            : nlopt::LN_BOBYQA;

    nlopt::opt solver(algorithm, x.size());

    solver.set_lower_bounds(lower);
    solver.set_upper_bounds(upper);
    solver.set_xtol_rel(options.xtol_rel);
    solver.set_ftol_rel(options.ftol_rel);
    solver.set_maxeval(
        static_cast<int>(options.max_evaluations)
    );

    solver.set_min_objective(
        &objective_trampoline,
        static_cast<void*>(&data)
    );

    double minimum = 0.0;

    nlopt::result nlopt_result = nlopt::FAILURE;

    try {
        nlopt_result = solver.optimize(x, minimum);
    } catch (const std::exception& error) {
        Result result;
        result.status = Status::Failure;
        result.converged = false;
        result.parameters = data.full;
        result.objective_value = minimum;
        result.message = error.what();
        return result;
    }

    Result result;
    result.parameters = data.full;
    result.objective_value = minimum;

    switch (nlopt_result) {
        case nlopt::SUCCESS:
        case nlopt::FTOL_REACHED:
        case nlopt::XTOL_REACHED:
        case nlopt::STOPVAL_REACHED:
            result.status = Status::Success;
            result.converged = true;
            result.message = "converged";
            break;

        case nlopt::MAXEVAL_REACHED:
        case nlopt::MAXTIME_REACHED:
            result.status = Status::MaxEvaluations;
            result.converged = false;
            result.message = "evaluation budget exhausted";
            break;

        case nlopt::INVALID_ARGS:
            result.status = Status::Invalid;
            result.converged = false;
            result.message = "invalid arguments";
            break;

        default:
            result.status = Status::Failure;
            result.converged = false;
            result.message = solver.get_errmsg();
            break;
    }

    result.evaluations = static_cast<std::size_t>(
        solver.get_numevals()
    );

    return result;
}

} // namespace optimization
} // namespace agnsf
