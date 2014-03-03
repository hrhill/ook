#ifndef OOK_CALL_SELECTOR_H_
#define OOK_CALL_SELECTOR_H_

#include <tuple>

namespace ook{
namespace detail{

// Meta function to select the right function call
// based on the properties of the return type.
template <typename F, typename X, typename State, int dim>
struct call_selector{};
/*
template <typename F, typename X, typename State>
struct call_selector<F, X, State, 1>{
    static
    void
    call(F f, const X& x, State& s){
        std::tie(s.fx) = f(x);
    }
};
*/
template <typename F, typename X, typename State>
struct call_selector<F, X, State, 2>{
    static
    void
    call(F f, const X& x, State& s){
        std::tie(s.fx, s.dfx) = f(x);
    }
};

template <typename F, typename X, typename State>
struct call_selector<F, X, State, 3>{
    static
    void
    call(F f, const X& x, State& s){
        std::tie(s.fx, s.dfx, s.H) =  f(x);
    }
};

} // ns detail

} // ns ook

#endif
