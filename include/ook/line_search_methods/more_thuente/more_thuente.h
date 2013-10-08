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

//    std::cout << "Entering more_thuente with arguments\n" << phi0 << ", " << dphi0 << ", " << a << std::endl;
    T phia, dphia;
    std::tie(phia, dphia) = phi(a);

    more_thuente_searcher<T> search(phi0, dphi0, a, opts.stpmax - opts.stpmin);

    message msg = message::start;
    uint attempts = 0;
    do{
        if(strong_wolfe_conditions(phi0, phia, dphi0, dphia, a, opts.ftol, opts.gtol))
        {
            msg = message::convergence;
            break;
        }
        std::tie(msg, a) = search(a, phia, dphia, opts);
        ++attempts;
/*
        if (attempts == opts.max_line_search_attempts){
            value = message::warning_max_line_search_attempts_reached;
            break;
        }
*/
        if (msg != message::update)
            break;

        std::tie(phia, dphia) = phi(a);
    }while (true);

    return std::make_tuple(msg, a);
}

} // ns line_search

} // ns ook

#endif