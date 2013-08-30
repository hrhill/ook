#ifndef OOK_OPTIMISE_H_
#define OOK_OPTIMISE_H_

#include "./schemes/state.h"

namespace ook{

/// main optimisation loop
template <typename Scheme, typename State, typename F, typename X>
State
optimise(F objective_function, const X& x0, const typename Scheme::options_type&
 options)
{
    //State s0 = Scheme::initialise(objective_function, x0, State());
    // Evaluate at initial point
    State s0;
    s0.xk = x0;
    auto fx_dfx = objective_function(x0);
    s0.fxk = std::get<0>(fx_dfx);
    s0.dfxk = std::get<1>(fx_dfx);

    do {
        // Choose descent direction
        s0.pk = -s0.dfxk;
        s0.alpha = 1.0;

        // do line search
        //State sk = Scheme::iterate(objective_function, s0);
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
