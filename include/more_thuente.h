#ifndef OOK_MORE_THUENTE_H_
#define OOK_MORE_THUENTE_H_

#include <tuple>
#include <iomanip>

#include "line_search_conditions.h"
#include "state.h"
#include "dcsrch.h"

namespace ook{

template <typename F>
state
more_thuente_line_search(F phi, state sk, const options& opts){

    std::tie(sk.fxap, sk.dfxap_dot_p) = phi(sk.a);
    sk.value = state_value::start;
    dcsrch_struct dcsrch_(sk.fx, sk.dfx_dot_p, sk.a, opts.stpmax - opts.stpmin);

    do{
        if (strong_wolfe_conditions(sk.fx, sk.fxap, sk.dfx_dot_p, sk.dfxap_dot_p, sk.a, opts.ftol, opts.gtol))
        {
            sk.value = state_value::convergence;
            break;
        }
        std::tie(sk.value, sk.a) = dcsrch_(sk.a, sk.fxap, sk.dfxap_dot_p, opts);
        if (sk.value != state_value::update)
            break;

        std::tie(sk.fxap, sk.dfxap_dot_p) = phi(sk.a);
    }while (true);

    return sk;
}

} // ns ook

#endif