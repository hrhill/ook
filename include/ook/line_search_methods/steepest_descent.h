#ifndef OOK_LINE_SEARCH_METHODS_STEEPEST_DESCENT_H_
#define OOK_LINE_SEARCH_METHODS_STEEPEST_DESCENT_H_

#include <tuple>
#include "ook/line_search_methods/line_search_method.h"

namespace ook{

namespace detail{

template <typename X>
struct steepest_descent{
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
        return s;
    }
    
    static
    vector_type
    descent_direction(state_type& s)
    {
        return -s.dfx;
    }

    static
    state_type
    update(state_type s){
        return s;
    }    
};

} // ns detail

template <typename F, typename X, typename Options>
std::tuple<ook::state_value, X>
steepest_descent(F objective_function, const X& x0, const Options& opts)
{
    return line_search_method<detail::steepest_descent<X>>(objective_function, x0, opts);
}

} //ns ook


#endif
