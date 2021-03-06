// Copyright 2013 Harry Hill
//
// This file is part of ook.
//
// ook is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ook is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with ook.  If not, see <http://www.gnu.org/licenses/>.

#ifndef OOK_FINITE_DIFFERENCES_BACKWARD_DIFFERENCE_HPP_
#define OOK_FINITE_DIFFERENCES_BACKWARD_DIFFERENCE_HPP_

#include <cmath>
#include <limits>
#include <tuple>

#include <algorithm>

#include "ook/finite_differences/detail/transform.hpp"

namespace ook
{
namespace finite_differences
{

/// \brief Backward difference approximation.
struct backward_difference
{
    /// \brief Calculate the backward difference approximation to the gradient
    /// of f.
    template <typename F, typename X>
    static auto
    gradient(F f, X x);

    /// \brief Calculate the backward difference approximation to the hessian of
    /// f.
    template <typename M, typename F, typename X>
    static auto
    hessian(F f, const X& x);
};

template <typename F, typename X>
auto
backward_difference::gradient(F f, X x)
{
    // Generate a set of sample points
    // (x + he_1, x + he_2, ..., x)
    const size_t n = std::distance(x.begin(), x.end());
    const double hmin = sqrt(std::numeric_limits<double>::epsilon());
    std::vector<X> sample_points(n + 1);
    X h(n);

    for (size_t i = 0; i < n; ++i)
    {
        const double xi = x[i];
        const double hx = hmin * (1 + fabs(xi));
        h[i] = hx;
        x[i] = xi - hx;
        sample_points[i] = x;
        x[i] = xi;
    }
    sample_points[n] = x;

    // evaluate function at each point
    std::vector<double> function_values(sample_points.size());
    ook::finite_differences::detail::transform(
        sample_points, function_values, f);

    // assemble
    X df(n);
    const double fx = function_values[n];
    for (size_t i = 0; i < n; ++i)
    {
        df[i] = (fx - function_values[i]) / h[i];
    }
    return std::make_tuple(fx, df);
}

template <typename M, typename F, typename X>
auto
backward_difference::hessian(F f, const X& x)
{
    const double eps = std::numeric_limits<double>::epsilon();
    const double hmin = std::pow(eps, 1.0 / 3.0);

    const size_t n = std::distance(x.begin(), x.end());

    X xi(x);
    X xj(x);

    const double fx = f(x);
    M H(n, n);

    for (size_t i = 0; i < n; ++i)
    {
        const double xii = xi[i];
        const double hi = hmin * (1 + fabs(xii));

        // f(x - h_i e_i)
        xi[i] -= hi;
        const double fxi = f(xi);
        xi[i] = xii;

        for (size_t j = 0; j <= i; ++j)
        {
            const double xjj = xj[j];
            const double hj = hmin * (1 + fabs(xjj));

            // f(x - h_j e_j)
            xj[j] -= hj;
            const double fxj = f(xj);

            // f(x - h_i e_i - h_j e_j)
            xj[i] -= hi;
            const double fxij = f(xj);

            H(i, j) = H(j, i) = (fxij - fxi - fxj + fx) / (hi * hj);

            xj[i] = xii;
            xj[j] = xjj;
        }
    }
    return std::make_tuple(fx, H);
}

} // ns finite_differences
} // ns ook

#endif
