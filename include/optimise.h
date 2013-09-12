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
        if (Scheme::termination_criterion(sk, options))
            break;
        
    } while(false);

    return sk;
}

} //ns ook

#endif
