#ifndef OOK_LINE_SEARCH_METHODS_FLETCHER_REEVES_H_
#define OOK_LINE_SEARCH_METHODS_FLETCHER_REEVES_H_

#include <tuple>
#include <algorithm>

#include "line_search_method.h"

namespace ook{

namespace detail{

template <typename X>
struct fletcher_reeves{
    typedef X vector_type;
    typedef typename X::value_type value_type;
    typedef state<X> state_type;    

    template <typename F>
    static
    state_type
    initialise(F objective_function, const X& x0)
    {
        state_type s(x0.size());
        std::tie(s.fx, s.dfx) = objective_function(x0);
        s.dfx0 = s.dfx;
        return s;
    }
    
    static
    vector_type
    descent_direction(state_type& s)
    {
        s.p = -s.dfx + s.beta * s.p;
        return s.p;
    }

    static
    state_type
    update(state_type s)
    {
        s.beta = detail::inner_product(s.dfx, s.dfx)/detail::inner_product(s.dfx0, s.dfx0);
        s.dfx0 = s.dfx;
        return s;
    }
};

} // ns detail

template <typename F, typename X, typename Options>
std::tuple<ook::state_value, X>
fletcher_reeves(F objective_function, const X& x0, const Options& opts)
{
    return line_search_method<detail::fletcher_reeves<X>>(objective_function, x0, opts);
}

} //ns ook


#endif