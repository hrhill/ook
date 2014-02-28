#ifndef OOK_OPTIONS_H_
#define OOK_OPTIONS_H_

#include <cassert>
#include <limits>

namespace ook{

template <typename T>
struct options
{
    options()
    :
        ftol(T(1e-03)),
        gtol(T(1e-01)),
        xtol(std::numeric_limits<T>::epsilon()),
        stpmin(std::numeric_limits<T>::epsilon()),
        stpmax(T(5.0)),
        max_function_evaluations(2000)
    {}

    options(const T& ftol_, const T& gtol_, const T& xtol_,
            const T& stpmin_, const T& stpmax_)
    :
        ftol(ftol_), gtol(gtol_), xtol(xtol_), stpmin(stpmin_), stpmax(stpmax_)
    {
        assert(ftol > T(0) && "ftol <= 0.");
        assert(gtol > T(0) && "gtol <= 0.");
        assert(xtol > T(0) && "xtol <= 0.");
        assert(stpmin >= T(0.0) && "stpmin < 0.");
        assert(stpmax > stpmin && "stpmax <= stpmin.");
    }

    // tolerance for the sufficient decrease condition.
    T ftol;
    // tolerance for the curvature condition.
    T gtol;
    // relative tolerance for an acceptable step.
    T xtol;
    // lower bound for the step.
    T stpmin;
    // upper bound for the step.
    T stpmax;

    uint max_function_evaluations;
};

} // ns ook

#endif
