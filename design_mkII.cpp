// User inputs
// - Objective function
// - Initial point to start iteration
// - Convergence criteria options
// 
// Outputs
// - Solution to optimisation problem
// - Gradient check functionality
// - Report on solver progress and final solution
#include <tuple>

struct real_type{};
struct vector_type{};
struct options_type{};

enum class optimiser_status{
    keep_going,
    problem_with_options,
    problem_with_objective_function,
    max_iterations_reached,
    no_step_possible
};

optimiser_status
validate_options(const options_type& options)
{
    return optimiser_status::keep_going;
}

template <typename ObjectiveFunction>
optimiser_status
validate_objective_function(ObjectiveFunction objective_function)
{
    return optimiser_status::keep_going;
}

template <typename Scheme, typename ObjectiveFunction, typename Domain, typename Observer>
std::tuple<optimiser_status, real_type, vector_type>
optimise(ObjectiveFunction objective_function, Domain x0, options_type& options, Observer observer)
{
    // Check options,]. restore them to sensible values, write them to the observer
    if (validate_options(options) != optimiser_status::keep_going)
        break

    // Validate the objective function if required
    if (validate_objective_function(objective_function) != optimiser_status::keep_going)
        break;

    // Setup main loop
    unsigned int iteration = 1;
    auto fx = objective_function(x0);

    do{
        // Choose descent direction
        auto p = Scheme::search_direction(fx);

        // Perform line search to get new point
        auto fxap = Scheme::line_search(fx, p);

        // Check for convergence at new point
        bool converged = true;
        if (converged)
            break;

        // Update
        fx = fxap;
        ++iteration;
    }while(true);

    return std::make_tuple(status, fxk, xk);
}
