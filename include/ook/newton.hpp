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

#include "linalg.hpp"

#include "ook/line_search_method.hpp"

namespace ook{
namespace detail{

/// \brief Take a matrix in LD format and convert
/// it to a lower cholesky matrix.
template <typename Matrix>
Matrix
convert_to_cholesky(const Matrix& LD)
{
    int n = linalg::num_rows(LD);
    Matrix L(LD);
    for (int j = 0; j < n; ++j){
        const double di = sqrt(L(j, j));
        L(j, j) = di;
        for (int i = j + 1; i < n; ++i){
            L(i, j) *= di;
        }
    }
    return L;
}

/// \brief Solve the system Ax = b where A is a
/// symmetric positive definite matrix.
template <typename Matrix, typename Vector>
Vector
solve(Matrix A, const Vector& b)
{
    Matrix LD = linalg::factorisations::gmw81(A);
    Matrix L = convert_to_cholesky(LD);

    Matrix b1(linalg::size(b), 1);

    linalg::column(b1, 0) = b;
    linalg::potrs(L, b1);

    return linalg::column(b1, 0);
}

} // ns detail

/// \brief Implementation of the required steps of line_search_method
/// for Newtons method.
template <typename X>
struct newton_impl
{
    template <typename T>
    newton_impl(const T&){}

    /// \brief The descent direction for the Newton method
    /// is determined by find the solution p to the system
    /// \f[
    ///          H p = -\nabla f(x).
    /// \f]
    ///
    template <typename State>
    X
    descent_direction(const State& s)
    {
        return -detail::solve(s.H, s.dfx);
    }

    template <typename T>
    void
    update(const T&){}
};

/// \brief The Newton algorithm.
/// \details Implementation of the Newton algorithm using the generic line
/// search function.
template <typename F, typename X, typename Options, typename Observer>
std::tuple<ook::message, X>
newton(F f, const X& x0, const Options& opts, Observer& observer)
{
    typedef newton_impl<X> scheme;
    line_search_method<scheme> method;
    return method(f, x0, opts, observer);

}

} //ns ook

#endif
