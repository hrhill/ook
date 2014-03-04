#ifndef OOK_LINE_SEARCH_BACKTRACKING_H_
#define OOK_LINE_SEARCH_BACKTRACKING_H_

#include <tuple>
#include <exception>
#include <cassert>

#include "./conditions.h"

namespace ook{
namespace line_search{

struct backtracking{

    template <typename F, typename T, typename Options>
    static
    std::tuple<T, T, T>
    search(F phi, const T& phi0, const T& dphi0, T a, const Options& opts)
    {
        T phia, dphia;
        T rho(0.9); // Need to pass this as an option.

        while(true){
            if (fabs(a) <= std::numeric_limits<T>::epsilon())
                break;
            std::tie(phia, dphia) = phi(a);
            if (sufficient_decrease_condition(phia, phi0, opts.ftol, a, dphi0))
                break;
            a *= rho;
        }
        return std::make_tuple(a, phia, dphia);
    }
};


} // ns line_search
} // ns ook

#endif
