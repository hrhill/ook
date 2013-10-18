#ifndef OOK_LINE_SEARCH_MORE_THUENTE_H_
#define OOK_LINE_SEARCH_MORE_THUENTE_H_

#include <tuple>

#include "ook/line_search_methods/message.h"

#include "line_search_conditions.h"
#include "more_thuente_searcher.h"

namespace ook{
namespace line_search{

template <typename T>
T
safeguarded(T at, T al, T delta, T amax){
    return std::min(at + delta * (at - al), amax);
}

template <typename T>
std::tuple<T, T>
update_interval(T phi0, T dphi0, T at, T phiat, T dphiat, T al, T phial, T dphial, T au, T mu)
{
    const double psiat = phiat - phi0 - mu * dphi0 * at;
    const double psial = phial - phi0 - mu * dphi0 * al;
    if (psiat > psial){
        return std::make_pair(al, at);
    }else{
        if (dphiat * (al - at) > 0){
            // Choose new at
            std::make_pair(at, au);
        }else{
            return std::make_pair(at, al);
        }
    }
}
/*
template <typename F, typename T, typename Options>
std::tuple<message, T>
more_thuente_v2(F phi, T phi0, T dphi0, T a, const Options& opts)
{
    message msg = message::start;

    if (dphi0 >= 0){
        msg = message::search_direction_is_not_a_descent_direction;
        return std::make_tuple(msg, a);
    }

    T phia, dphia, phial, dphial;
    do{
        std::tie(phia, dphia) = phi(a);        

        if(strong_wolfe_conditions(phi0, phia, dphi0, dphia, a, opts.ftol, opts.gtol))
        {
            msg = message::convergence;
            break;
        }
        // Update interval
        std::tie(al, au) = update_interval(phi0, dphi0, a, phia, dphia, al, phial, dphial, au, mu);
        // Choose safeguarded interval
        a = safeguarded(a, al, delta, amax);
    }while (true);

    return std::make_tuple(msg, a);
}
*/

template <typename F, typename T, typename Options>
std::tuple<message, T>
more_thuente(F phi, T phi0, T dphi0, T a, const Options& opts)
{
    message msg = message::start;

    if (dphi0 >= 0){
        msg = message::search_direction_is_not_a_descent_direction;
        return std::make_tuple(msg, a);
    }

    T phia, dphia;
    std::tie(phia, dphia) = phi(a);

    more_thuente_searcher<T, Options> search(phi0, dphi0, a, opts.stpmax - opts.stpmin, opts);

    do{
        if (isnan(a)){
            a = opts.stpmax;
            msg = message::error_step_greater_than_stpmax;
            break;
        }
        if(strong_wolfe_conditions(phi0, phia, dphi0, dphia, a, opts.ftol, opts.gtol))
        {
            msg = message::convergence;
            break;
        }
        std::tie(msg, a) = search(a, phia, dphia);
        if (msg != message::update)
            break;

        if (a <= opts.stpmin){
            a = opts.stpmin;
            msg = message::error_step_less_than_stpmin;
            break;
        }
        if (a >= opts.stpmax){
            a = opts.stpmax;
            msg = message::error_step_greater_than_stpmax;
            break;            
        }

        std::tie(phia, dphia) = phi(a);
    }while (true);

    return std::make_tuple(msg, a);
}

} // ns line_search

} // ns ook

#endif