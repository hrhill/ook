#ifndef OOK_OPTIMISE_H_
#define OOK_OPTIMISE_H_

#include "./schemes/state.h"

namespace ook{

/// main optimisation loop
template <typename Scheme, typename State, typename F, typename X>
State
optimise(F objective_function, const X& x0, const typename Scheme::options_type& options)
{
    State sk = Scheme::initialise(objective_function, x0, State());
    // Evaluate at initial point
    sk.x = x0;
    auto fx_dfx = objective_function(x0);
    sk.fx = std::get<0>(fx_dfx);
    sk.dfx = std::get<1>(fx_dfx);

    do {
        // Choose descent direction
        vector_type p = Scheme::get_descent_direction(objective_function, sk);
        // do line search
        sk = Scheme::search_along_direction(objective_function, sk);
        //s0 = Scheme::check_and_advance(s0, sk, options);

    } while(false);

    return s0;
}

/*
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

*/
}

#endif
