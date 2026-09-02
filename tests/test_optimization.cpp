#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include <optimization/optimization.hpp>

namespace {

bool close(
    double a,
    double b,
    double rtol = 1e-4,
    double atol = 1e-4
)
{
    return std::abs(a - b) <= atol + rtol * std::abs(b);
}


class Quadratic : public agnsf::optimization::Objective {
public:
    Quadratic(double a, double b)
        : a_(a), b_(b)
    {
    }

    double evaluate(
        const std::vector<double>& x
    ) const override
    {
        const double dx = x[0] - a_;
        const double dy = x[1] - b_;
        return dx * dx + dy * dy;
    }

private:
    double a_;
    double b_;
};


void test_recover_quadratic()
{
    agnsf::optimization::Options options;

    std::vector<agnsf::optimization::Parameter> parameters = {
        {0.0, -10.0, 10.0},
        {0.0, -10.0, 10.0},
    };

    const Quadratic objective(2.0, -3.0);

    const agnsf::optimization::Result result =
        agnsf::optimization::minimize(objective, parameters, options);

    assert(result.converged);
    assert(result.status == agnsf::optimization::Status::Success);
    assert(close(result.parameters[0], 2.0));
    assert(close(result.parameters[1], -3.0));
    assert(result.objective_value < 1e-6);
}


void test_solution_at_bound()
{
    agnsf::optimization::Options options;

    // Minimum of (x - 5)^2 lies outside [0, 3]; the solution must be
    // pinned at the active bound.
    std::vector<agnsf::optimization::Parameter> parameters = {
        {0.0, 0.0, 3.0},
    };

    const Quadratic objective(5.0, 0.0);

    const agnsf::optimization::Result result =
        agnsf::optimization::minimize(objective, parameters, options);

    assert(result.converged);
    assert(close(result.parameters[0], 3.0, 1e-6, 1e-6));
    assert(close(result.objective_value, 4.0, 1e-6, 1e-6));
}


void test_fixed_parameter()
{
    agnsf::optimization::Options options;

    // y is pinned at 2.0; only x is optimized.
    std::vector<agnsf::optimization::Parameter> parameters = {
        {0.0, -10.0, 10.0},
        {2.0, -10.0, 10.0, true},
    };

    const Quadratic objective(1.0, 5.0);

    const agnsf::optimization::Result result =
        agnsf::optimization::minimize(objective, parameters, options);

    assert(result.converged);
    assert(close(result.parameters[0], 1.0));
    assert(result.parameters[1] == 2.0);
    assert(close(result.objective_value, 9.0, 1e-6, 1e-6));
}


void test_degenerate_bounds_pin()
{
    agnsf::optimization::Options options;

    // lower == upper pins the parameter without the fixed flag.
    std::vector<agnsf::optimization::Parameter> parameters = {
        {0.0, -10.0, 10.0},
        {2.0, 2.0, 2.0},
    };

    const Quadratic objective(1.0, 5.0);

    const agnsf::optimization::Result result =
        agnsf::optimization::minimize(objective, parameters, options);

    assert(result.converged);
    assert(close(result.parameters[0], 1.0));
    assert(result.parameters[1] == 2.0);
}


void test_all_fixed()
{
    agnsf::optimization::Options options;

    std::vector<agnsf::optimization::Parameter> parameters = {
        {1.0, 0.0, 2.0, true},
        {-2.0, -5.0, 0.0, true},
    };

    const Quadratic objective(0.0, 0.0);

    const agnsf::optimization::Result result =
        agnsf::optimization::minimize(objective, parameters, options);

    assert(result.converged);
    assert(result.evaluations == 1);
    assert(result.parameters[0] == 1.0);
    assert(result.parameters[1] == -2.0);
    assert(close(result.objective_value, 5.0));
}


void test_initial_value_clamped()
{
    agnsf::optimization::Options options;

    // Initial value outside the box is clamped before optimizing.
    std::vector<agnsf::optimization::Parameter> parameters = {
        {50.0, -10.0, 10.0},
        {-50.0, -10.0, 10.0},
    };

    const Quadratic objective(2.0, -3.0);

    const agnsf::optimization::Result result =
        agnsf::optimization::minimize(objective, parameters, options);

    assert(result.converged);
    assert(close(result.parameters[0], 2.0));
    assert(close(result.parameters[1], -3.0));
}


void test_invalid_bounds_throw()
{
    agnsf::optimization::Options options;

    std::vector<agnsf::optimization::Parameter> parameters = {
        {0.0, 5.0, -5.0}, // lower > upper
    };

    const Quadratic objective(0.0, 0.0);

    bool threw = false;

    try {
        agnsf::optimization::minimize(objective, parameters, options);
    } catch (const std::invalid_argument&) {
        threw = true;
    }

    assert(threw);

    // Non-finite bounds.
    parameters[0] = {0.0, -std::numeric_limits<double>::infinity(), 1.0};
    threw = false;

    try {
        agnsf::optimization::minimize(objective, parameters, options);
    } catch (const std::invalid_argument&) {
        threw = true;
    }

    assert(threw);

    // Empty parameter list.
    threw = false;

    try {
        agnsf::optimization::minimize(objective, {}, options);
    } catch (const std::invalid_argument&) {
        threw = true;
    }

    assert(threw);
}


void test_evaluation_budget()
{
    agnsf::optimization::Options options;
    options.algorithm = agnsf::optimization::Algorithm::NelderMead;
    options.max_evaluations = 4;
    options.xtol_rel = 1e-12;
    options.ftol_rel = 1e-12;

    std::vector<agnsf::optimization::Parameter> parameters = {
        {0.0, -10.0, 10.0},
        {0.0, -10.0, 10.0},
    };

    const Quadratic objective(2.0, -3.0);

    const agnsf::optimization::Result result =
        agnsf::optimization::minimize(objective, parameters, options);

    assert(result.status ==
        agnsf::optimization::Status::MaxEvaluations);
    assert(!result.converged);
    assert(result.evaluations <= options.max_evaluations);
}

} // namespace


int main()
{
    test_recover_quadratic();
    test_solution_at_bound();
    test_fixed_parameter();
    test_degenerate_bounds_pin();
    test_all_fixed();
    test_initial_value_clamped();
    test_invalid_bounds_throw();
    test_evaluation_budget();

    return 0;
}
