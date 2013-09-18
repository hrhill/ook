#ifndef OOK_MORE_THUENTE_H_
#define OOK_MORE_THUENTE_H_

#include <tuple>

#include "line_search_conditions.h"
#include "state_value.h"
#include "dcsrch.h"

namespace ook{

template <typename F, typename T, typename Options>
std::tuple<state_value, T>
more_thuente_line_search(F phi, T phi0, T dphi0, T a, const Options& opts){

    T phia, dphia;
    std::tie(phia, dphia) = phi(a);
    dcsrch_struct<T> dcsrch_(phi0, dphi0, a, opts.stpmax - opts.stpmin);

    state_value value = state_value::start;

    do{
        if(strong_wolfe_conditions(phi0, phia, dphi0, dphia, a, opts.ftol, opts.gtol))
        {
            value = state_value::convergence;
            break;
        }
        std::tie(value, a) = dcsrch_(a, phia, dphia, opts);
        if (value != state_value::update)
            break;

        std::tie(phia, dphia) = phi(a);
    }while (true);

    return std::make_tuple(value, a);
}

} // ns ook

#endif