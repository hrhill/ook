#ifndef OPTIM_SCHEMES_BFGS_H_
#define OPTIM_SCHEMES_BFGS_H_

namespace ook{

struct newton_options{};

struct newton{

    typedef newton_options options_type;

    template <typename F, typename X, typename State>
    static
    State
    update(F f, const X& x)
    {
        return State();
    }

    template <typename F, typename X, typename State>
    static
    State
    initialise(F f, const X& x, const State&)
    {
        return update<F, X, State>(f, x);        
    }

    template <typename F, typename State>
    static 
    State
    iterate(F f, const State& s0)
    {
        return s0;
    }

    template <typename State>
    static
    State
    check_and_advance(const State& s_old, const State& s_new, const options_type& opts)
    {
        // Check relative change in x
        // Check relative change in fx
        // Check first order conditions
        // If we haven't converged, move s_new into s_old
        // and return the 
        return s_old;
    }

};

} // hs optim

#endif
