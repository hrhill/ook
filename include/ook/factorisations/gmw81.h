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

#ifndef OOK_FACTORISATIONS_GMW81_H_
#define OOK_FACTORISATIONS_GMW81_H_

#include <cassert>
#include <algorithm>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "ook/factorisations/tools.h"

namespace ook{

namespace factorisations{

namespace detail{

template <typename Matrix>
typename Matrix::value_type
calculate_theta(const Matrix& c, const int j){

    const int n = c.size1();
    auto jth_col = boost::numeric::ublas::matrix_column<const Matrix>(c, j);
    typename Matrix::value_type theta(0.0);
    if (j < n )
        for (int i = j + 1; i < n; ++i){
            theta = std::max(theta, jth_col(i));
        }
    return theta;
}

} // ns detail

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
\tparam Matrix Templated by matrix type, but currently restricted to Boost.Ublas.
**/
template <typename Matrix>
Matrix
gmw81(Matrix G)
{
    using namespace boost::numeric::ublas;

    typedef typename Matrix::value_type real_type;
    typedef typename Matrix::size_type size_type;

    // MC1
    const size_type n = G.size1();
    const real_type eps = std::numeric_limits<real_type>::epsilon();
    const real_type nu = std::max(1.0, sqrt(std::pow(n, 2) - 1.0));
    const real_type gamma = std::get<1>(tools::max_magnitude_diagonal(G));
    const real_type eta = tools::max_magnitude_off_diagonal(G);
    const real_type beta2 = std::max({gamma, eta/nu, eps});
    const real_type delta = sqrt(eps);

    // MC2
    Matrix c(n, n, 0);
    Matrix L(n, n, 0);

    for (size_type i = 0; i < n; ++i){
        c(i, i) = G(i, i);
    }
    for (size_type j = 0; j < n; ++j){
        // MC3
        real_type cqq;
        size_type q;
        std::tie(q, cqq) = tools::max_magnitude_diagonal(matrix_range<Matrix>(c, range(j, n), range(j, n)));

        // Pivoting
        if (q != j){
            matrix_row<Matrix> rowq(G, q);
            matrix_row<Matrix> rowj(G, j);
            rowq.swap(rowj);

            matrix_row<Matrix> colq(G, q);
            matrix_row<Matrix> colj(G, j);
            colq.swap(colj);
        }

        // MC4
        for (size_type s = 0; s < j; ++s){
            L(j, s) = c(j, s)/L(s, s);
        }
        for (size_type i = j + 1; i < n; ++i){
            c(i, j) = G(i, j);
            for (int s = 0; s < j; ++s){
                c(i, j) -= L(j, s) * c(i, s);
            }
        }
        real_type theta = detail::calculate_theta(c, j);
        // MC 5
        L(j, j) = std::max({delta, fabs(c(j, j)), std::pow(theta, 2)/beta2});

        for (size_type i = j + 1; i < n; ++i){
            c(i, i) -= std::pow(c(i, j), 2)/L(j, j);
        }
    }
    return L;
}

} // ns factorisation

} // ns ook

#endif
