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

#include "ook/state.hpp"
#include "ook/line_search/more_thuente.hpp"
#include "ook/line_search_method.hpp"

namespace ook{
namespace detail{

template <typename X>
struct newton_state
{
    typedef X vector_type;
    typedef typename std::remove_reference<decltype(X()[0])>::type value_type;
    typedef typename linalg::associated_matrix<X>::type matrix_type;

    newton_state(const int n = 0)
    :
        fx(0),
        dfx(n),
        dfx0(n),
        p(n),
        dx(n),
        H(n, n, 0.0),
        a(1),
        iteration(0),
        nfev(0),
        tag(state_tag::init)
    {
        for (int i = 0; i < n; ++i){
            H(i, i) = 1.0;
        }
    }

    value_type fx;
    vector_type dfx;
    vector_type dfx0;
    vector_type p;
    vector_type dx;
    matrix_type H;
    value_type a;
    int iteration;
    int nfev;
    state_tag tag;
    message msg;
};

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

/// \brief Implementation of the required steps of line_search_method
/// for Newtons method.
template <typename X>
struct newton{
    typedef X vector_type;
    typedef typename X::value_type value_type;
    typedef newton_state<X> state_type;

    template <typename F>
    state_type
    initialise(F objective_function, const X& x0)
    {
        state_type s(linalg::size(x0));
        std::tie(s.fx, s.dfx, s.H) = objective_function(x0);
        return s;
    }

    state_type
    descent_direction(state_type s)
    {
        ++s.iteration;
        s.p = -detail::solve(s.H, s.dfx);
        return s;
    }

    state_type
    update(const state_type& s)
    {
        return s;
    }

    ook::line_search::more_thuente search;
};

} // ns detail

/// \brief The Newton algorithm.
/// \details Implementation of the Newton algorithm using the generic line
/// search function.
template <typename F, typename X, typename Options, typename Observer>
std::tuple<ook::message, X>
newton(F f, const X& x0, const Options& opts, Observer& observer)
{
    typedef detail::newton<X> scheme;
    line_search_method<scheme> method;
    return method(f, x0, opts, observer);

}

} //ns ook

#endif
