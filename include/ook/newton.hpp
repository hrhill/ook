
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

#ifndef OOK_NEWTON_HPP_
#define OOK_NEWTON_HPP_

#include <tuple>

#include "ook/line_search/mcsrch.hpp"
#include "ook/line_search_method.hpp"

namespace ook
{
namespace detail
{

double
calculate_theta(const matrix& c, const size_t j)
{
    const size_t n = c.rows();
    double theta(0.0);
    if (j < n)
    {
        for (size_t i = j + 1; i < n; ++i)
        {
            theta = std::max(theta, c(i, j));
        }
    }
    return theta;
}

/// \brief Get the value of maximum magnitude on the diagonal of the matrix.
std::tuple<int, double>
max_magnitude_diagonal(const matrix& m)
{
    const int n = m.rows();

    if (n == 0)
        return std::make_tuple(0, 0.0);

    int idx = 0;
    double mmd = fabs(m(0, 0));
    for (int i = 1; i < n; ++i)
    {
        const double mii = fabs(m(i, i));
        idx = mmd < mii ? i : idx;
        mmd = idx == i ? mii : mmd;
    }
    return std::make_tuple(idx, mmd);
}

/// \brief Get the value of maximum magnitude off the diagonal of the matrix.
double
max_magnitude_off_diagonal(const matrix& m)
{
    const int n = m.rows();

    double mmd = fabs(m(0, 0));
    for (int j = 0; j < n; ++j)
    {
        for (int i = j + 1; i < n; ++i)
        {
            mmd = std::max(fabs(m(i, j)), mmd);
        }
    }
    return mmd;
}

/// \brief Modified cholesky factorisation of Gill, Murray and Wright.
/**
\details Modified cholesky factorisation of Gill, Murray and Wright
from Ch.4 page 111 of "Practical Optimisation".
The algorithm is commonly referred to as gmw81, hence the name here.
The function takes a symmetric matrix G, and returns a modified variant
of the LDLt factorisation. The returned matrix has layout,
\f[
\begin{pmatrix}
    d_{11} & 0      & \cdots  & 0 \\
    l_{21} & d_{22} & \cdots  & \vdots\\
    \vdots &        & \ddots  & \vdots\\
    l_{nn} & l_{n2} & \cdots   & d_{nn}\\
    \end{pmatrix}
\f]
**/
void
gmw81(matrix& G)
{
    // MC1
    const size_t n = G.rows();
    const double eps = std::numeric_limits<double>::epsilon();
    const double nu = std::max(1.0, sqrt(std::pow(n, 2) - 1.0));
    const double gamma = std::get<1>(max_magnitude_diagonal(G));
    const double eta = max_magnitude_off_diagonal(G);
    const double beta2 = std::max({gamma, eta / nu, eps});
    const double delta = sqrt(eps);

    // MC2
    matrix c(n, n, 0);

    for (size_t i = 0; i < n; ++i)
    {
        c(i, i) = G(i, i);
    }
    for (size_t j = 0; j < n; ++j)
    {
        // MC3
        /*
                double cqq;
                size_t q;
                std::tie(q, cqq) =
           max_magnitude_diagonal(matrix_range<matrix>(c, range(j, n), range(j,
           n)));

                // Pivoting
                if (q != j){
                    matrix_row<matrix> rowq(G, q);
                    matrix_row<matrix> rowj(G, j);
                    rowq.swap(rowj);

                    matrix_column<matrix> colq(G, q);
                    matrix_column<matrix> colj(G, j);
                    colq.swap(colj);
                }
        */
        // MC4
        for (size_t s = 0; s < j; ++s)
        {
            G(j, s) = c(j, s) / G(s, s);
        }
        for (size_t i = j + 1; i < n; ++i)
        {
            c(i, j) = G(i, j);
            for (size_t s = 0; s < j; ++s)
            {
                c(i, j) -= G(j, s) * c(i, s);
            }
        }
        double theta = detail::calculate_theta(c, j);
        // MC 5
        G(j, j) = std::max({delta, fabs(c(j, j)), std::pow(theta, 2) / beta2});

        for (size_t i = j + 1; i < n; ++i)
        {
            c(i, i) -= std::pow(c(i, j), 2) / G(j, j);
        }
    }
}

/// \brief Take a matrix in LD format and convert
/// it to a lower cholesky matrix.
inline void
convert_to_cholesky(matrix& L)
{
    int n = L.rows();
    for (int j = 0; j < n; ++j)
    {
        const double di = sqrt(L(j, j));
        L(j, j) = di;
        for (int i = j + 1; i < n; ++i)
        {
            L(i, j) *= di;
        }
    }
}

/// \brief Solve the system Ax = b where A is a
/// symmetric positive definite matrix.
vector
solve(matrix A, vector b)
{
    gmw81(A);
    convert_to_cholesky(A);
    potrs(A, b, 'L');

    return b;
}

} // ns detail

/// \brief Implementation of the required steps of line_search_method
/// for Newtons method.
struct newton_impl
{
    struct state
    {
        matrix H;
    };

    template <typename T>
    explicit newton_impl(const T&)
    {
    }

    /// \brief The descent direction for the Newton method
    /// is determined by find the solution p to the system
    /// \f[
    ///          H p = -\nabla f(x).
    /// \f]
    ///
    template <typename State>
    vector
    descent_direction(const State& s)
    {
        return -detail::solve(s.H, s.dfx);
    }

    template <typename T>
    void
    update(const T&)
    {
    }
};

/// \brief The Newton algorithm.
/// \details Implementation of the Newton algorithm using the generic line
/// search function.
template <typename F, typename Options, typename Observer>
typename line_search_method<newton_impl, line_search::mcsrch>::state_type
newton(F f, const vector& x0, const Options& opts, Observer& observer)
{
    typedef newton_impl scheme;
    line_search_method<scheme, line_search::mcsrch> method;
    return method(f, x0, opts, observer);
}

} // ns ook

#endif
