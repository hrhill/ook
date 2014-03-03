#ifndef OOK_LINE_SEARCH_CONDITIONS_H_
#define OOK_LINE_SEARCH_CONDITIONS_H_

#include <cassert>
#include <cmath>

namespace ook{

namespace line_search{

template <typename T>
bool
sufficient_decrease_condition(const T& fxap,
                              const T& fx,
                              const T& c,
                              const T& a,
                              const T& dfx_dot_p)
{
    assert(T(0.0) < c && c < T(1.0));
    return fxap <= fx + c * a * dfx_dot_p;
}

template <typename T>
bool
curvature_condition(const T& dfxap_dot_p,
                    const T& c,
                    const T& dfx_dot_p)
{
    assert(T(0.0) < c && c < T(1.0));
    return fabs(dfxap_dot_p) <= c * fabs(dfx_dot_p);
}

template <typename T>
bool
strong_wolfe_conditions(const T& fx,
                        const T& fxap,
                        const T& dfx_dot_p,
                        const T& dfxap_dot_p,
                        const T& a,
                        const T& c1 = T(1e-04),
                        const T& c2 = T(9e-01))
{
    return sufficient_decrease_condition(fxap, fx, c1, a, dfx_dot_p)
            && curvature_condition(dfx_dot_p, dfxap_dot_p, c2);
}

} // ns line_search

} // ns ook

#endif
