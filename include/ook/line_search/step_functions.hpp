#ifndef OOK_LINE_SEARCH_STEP_FUNCTIONS_HPP_
#define OOK_LINE_SEARCH_STEP_FUNCTIONS_HPP_

#include <boost/math/special_functions/sign.hpp>

namespace ook{
namespace line_search{

/// \brief Return true if a is closer to c than b is.
/// \tparam T Numeric type, supporting operators +, *, / and -.
template <typename T>
bool
is_closer(T a, T b, T c)
{
    return fabs(a - c) < fabs(b - c);
}

/// \brief Calculate the secant step.
/// \tparam T Numeric type, supporting operators *, / and -.
/// \param x Point x.
/// \param dx Derivative of function at x.
/// \param y Point y.
/// \param dy Derivative of function at y.
template <typename T>
T
secant_step(T x, T dx, T y, T dy)
{
    return y - dy * (y - x) / (dy - dx);
}

/// \brief Calculate the value which minimizes the quadratic interpolant of
/// the input values.
/// \tparam T Numeric type, supporting operators +, *, / and -.
/// \param x Point x.
/// \param fx Function value at x.
/// \param dx Derivative of function at x.
/// \param y Point y.
/// \param fy Function value at y.
template <typename T>
T
quadratic_step(T x, T fx, T dx, T y, T fy)
{
    return x + dx / ((fx - fy) / (y - x) + dx) / 2.0 * (y - x);
}

/// \brief Calculate the value which minimizes the cubic interpolant of the
/// input values. The function normalizes values in the square root.
/// \tparam T Numeric type, supporting operators +, *, / and -.
/// \param x Point x.
/// \param fx Function value at x.
/// \param dx Derivative of function at x.
/// \param y Point y.
/// \param fy Function value at y.
/// \param dy Derivative of function at y.
template <typename T>
T
cubic_step(T x, T fx, T dx, T y, T fy, T dy)
{
    using boost::math::sign;

    const T d1 = dx + dy - 3 * (fx - fy)/(x - y);
    const T s = std::max({fabs(d1), fabs(dx), fabs(dy)});
    const T d2 = sign(y - x) * s * sqrt(pow(d1 / s, 2) - (dx / s) * (dy / s));
    return y - (y - x) * (dy + d2 - d1)/(dy - dx + 2.0 * d2);
}


} // ns line_search
} // ns ook

#endif