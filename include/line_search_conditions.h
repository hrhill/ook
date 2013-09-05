#ifndef OOK_LINE_SEARCH_CONDITIONS_H_
#define OOK_LINE_SEARCH_CONDITIONS_H_

#include <cassert>

namespace ook{

template <typename State>
bool
sufficient_decrease(const State& s, const typename State::real_type& c)
{
    return s.fxk_alpha_pk <= s.fxk + c * s.alpha * s.dfxk_dot_pk;    
}

template <typename State>
bool
curvature_condition(const State& s, const typename State::real_type& c)
{
    return s.dfxk_alpha_pk_dot_pk <= c * s.dfxk_dot_pk;    
}

template <typename State>
bool
auxillary_curvature_condition(const State& s, const typename State::real_type& c)
{
    return fabs(s.dfxk_alpha_pk_dot_pk) <= c * fabs(s.dfxk_dot_pk);    
}

template <typename State>
bool
step_length_control(const State& s, const typename State::real_type& c)
{
    return s.fxk + (1.0 - c) * s.alpha * s.dfxk_dot_pk;
}

template <typename State>
bool
wolfe_conditions(const State& s, 
                 const typename State::real_type& c1 = typename State::real_type(1e-04),
                 const typename State::real_type& c2 = typename State::real_type(9e-01))
{
    typedef typename State::real_type real_type;
    assert(real_type(0.0) < c1 && c1 < c2 && c2 < real_type(1.0));
    return sufficient_decrease(s, c1) && curvature_condition(s, c2);
}

template <typename State>
bool
strong_wolfe_conditions(const State& s, 
                 const typename State::real_type& c1 = typename State::real_type
(1e-04),
                 const typename State::real_type& c2 = typename State::real_type
(9e-01))
{
    typedef typename State::real_type real_type;
    assert(real_type(0.0) < c1 && c1 < c2 && c2 < real_type(1.0));
    return sufficient_decrease(s, c1) && auxilliary_curvature_condition(s, c2);
}

template<typename State>
bool
goldstein_conditions(const State& s, const typename State::real_type& c = typename State::real_type(0.25))
{
    typedef typename State::real_type real_type;
    assert(c > real_type(0.0) && c < real_type(1.0));
    return sufficient_decrease(s, c) && step_length_control(s, c);
}

}
#endif