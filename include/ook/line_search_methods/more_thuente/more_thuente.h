#ifndef OOK_LINE_SEARCH_MORE_THUENTE_H_
#define OOK_LINE_SEARCH_MORE_THUENTE_H_

#include <tuple>

#include "ook/line_search_methods/message.h"

#include "line_search_conditions.h"
#include "more_thuente_searcher.h"

namespace ook{
namespace line_search{

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
        if(strong_wolfe_conditions(phi0, phia, dphi0, dphia, a, opts.ftol, opts.gtol))
        {
            msg = message::convergence;
            break;
        }
        std::tie(msg, a) = search(a, phia, dphia);
        if (msg != message::update)
            break;

        std::tie(phia, dphia) = phi(a);
    }while (true);

    return std::make_tuple(msg, a);
}

} // ns line_search

} // ns ook

#endif