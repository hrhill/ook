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

#ifndef OOK_LINE_SEARCH_METHODS_BFGS_HPP_
#define OOK_LINE_SEARCH_METHODS_BFGS_HPP_

#include <tuple>

#include "linalg.hpp"

#include "ook/line_search_method.hpp"
#include "ook/line_search/more_thuente.hpp"

namespace ook{
namespace detail{

template <typename X>
struct
bfgs_state
{
    typedef X vector_type;
    typedef typename linalg::associated_matrix<X>::type matrix_type;

    bfgs_state(const X& dfx)
    :
        dfx(dfx),
        p(dfx.size()),
        H(dfx.size(), dfx.size(), 0.0)
    {
        const int n = dfx.size();
        for (int i = 0; i < n; ++i){
            H(i, i) = 1.0;
        }
    }

    vector_type dfx;
    vector_type p;
    matrix_type H;
};

/// \brief Implementation of the required steps of line_search_method
/// for BFGS method.
template <typename X>
struct bfgs
{
    typedef typename std::remove_reference<decltype(X()[0])>::type value_type;
    typedef typename linalg::associated_matrix<X>::type matrix_type;
    typedef bfgs_state<X> state_type;

    template <typename State>
    bfgs(const State& s)
    : state(s.dfx)
    {}

    template <typename State>
    X
    descent_direction(const State& s)
    {
        linalg::gemv(-1.0, state.H, s.dfx, 0.0, state.p);
        return state.p;
    }

    template <typename State>
    void
    update(const State& s)
    {
        X yk(s.dfx - state.dfx);
        const value_type rho = 1.0/linalg::inner_prod(yk, s.dx);

        const int n = linalg::size(s.dfx);
        matrix_type Z(linalg::identity_matrix<matrix_type>(n) - rho * linalg::outer_prod(s.dx, yk));

        matrix_type tmp(n, n);
        linalg::gemm(1.0, state.H, linalg::trans(Z), 0.0, tmp);
        linalg::gemm(1.0, Z, tmp, 0.0, state.H);
        matrix_type ss = rho * linalg::outer_prod(s.dx, s.dx);
        state.H += ss;
        state.dfx = s.dfx;
    }

    ook::line_search::more_thuente search;
    state_type state;
};

} // ns detail

/// \brief The Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm.
template <typename F, typename X, typename Options, typename Observer>
std::tuple<ook::message, X>
bfgs(F f, const X& x0, const Options& opts, Observer& observer)
{
    typedef detail::bfgs<X> scheme;
    line_search_method<scheme> method;
    return method(f, x0, opts, observer);
}

} //ns ook


#endif
