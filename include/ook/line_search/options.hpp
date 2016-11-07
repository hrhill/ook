#ifndef OOK_LINE_SEARCH_OPTIONS_HPP_
#define OOK_LINE_SEARCH_OPTIONS_HPP_

#include <limits>

namespace ook
{
namespace line_search
{

/// \brief Options structure for line searches.
template <typename T>
struct options
{
    options()
        : ftol(T(1e-03)),
          gtol(T(9e-01)),
          stpmin(std::numeric_limits<T>::epsilon()),
          stpmax(T(5.0)),
          maxfev(20)
    {
    }

    options(T ftol, T gtol, T stpmin, T stpmax)
        : ftol(ftol), gtol(gtol), stpmin(stpmin), stpmax(stpmax), maxfev(20)
    {
        assert(ftol > T(0) && "ftol <= 0.0");
        assert(gtol > T(0) && "gtol <= 0.0");
        assert(stpmin >= T(0.0) && "stpmin < 0.0");
        assert(stpmax > stpmin && "stpmax <= stpmin.");
    }

    /// \brief The constant used in the sufficient decrease condition.
    T ftol;

    /// \brief The constant used in the curvature condition.
    T gtol;

    /// \brief The minimum allowable step size.
    T stpmin;

    /// \brief The maximum step size
    T stpmax;

    /// \brief The maximum number of function evaluations allowed
    /// for a given line search.
    int maxfev;
};

} // ns line search
} // ns ook

#endif
