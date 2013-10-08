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
more_thuente(F phi, T phi0, T dphi0, T a, const Options& opts){

    T phia, dphia;
    std::tie(phia, dphia) = phi(a);
    more_thuente_searcher<T> search(phi0, dphi0, a, opts.stpmax - opts.stpmin);

    message value = message::start;
    uint attempts = 0;
    do{
        if(strong_wolfe_conditions(phi0, phia, dphi0, dphia, a, opts.ftol, opts.gtol))
        {
            value = message::convergence;
            break;
        }
        std::tie(value, a) = search(a, phia, dphia, opts);
        ++attempts;

        if (attempts == opts.max_line_search_attempts){
            a = T(0);
            value = message::warning_max_line_search_attempts_reached;
            break;
        }

        if (value != message::update)
            break;

        std::tie(phia, dphia) = phi(a);
    }while (true);

    return std::make_tuple(value, a);
}

} // ns line_search

} // ns ook

#endif