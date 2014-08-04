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

#ifndef OOK_LINE_SEARCH_METHODS_BFGS_H_
#define OOK_LINE_SEARCH_METHODS_BFGS_H_

#include <tuple>

#include "linalg.hpp"

#include "ook/state.h"
#include "ook/line_search_method.h"
#include "ook/line_search/more_thuente.h"

namespace ook{
namespace detail{

template <typename X>
struct
bfgs_state{
    typedef X vector_type;
    typedef typename std::remove_reference<decltype(X()[0])>::type value_type;

    typedef typename linalg::associated_matrix<X>::type matrix_type;

    bfgs_state(const int n = 0)
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

/// \brief Implementation of the required steps of line_search_method
/// for BFGS method.
template <typename X>
struct bfgs{
    typedef X vector_type;
    typedef typename std::remove_reference<decltype(X()[0])>::type value_type;
    typedef typename linalg::associated_matrix<X>::type matrix_type;
    typedef bfgs_state<X> state_type;

    template <typename F>
    state_type
    initialise(F objective_function, const X& x0)
    {
        state_type s(linalg::size(x0));
        std::tie(s.fx, s.dfx) = objective_function(x0);
        s.dfx0 = s.dfx;
        return s;
    }

    state_type
    descent_direction(state_type s)
    {
        ++s.iteration;
        linalg::gemv(-1.0, s.H, s.dfx, 0.0, s.p);
        return s;
    }

    state_type
    update(state_type s){
        X y(s.dfx - s.dfx0);
        const value_type rho = 1.0/linalg::inner_prod(y, s.dx);

        const int n = linalg::size(s.dfx);
        matrix_type Z(linalg::identity_matrix<matrix_type>(n) - rho * linalg::outer_prod(s.dx, y));
        matrix_type ss = rho * linalg::outer_prod(s.dx, s.dx);

        if (s.iteration == 1){
            const value_type hii = linalg::inner_prod(s.dx, s.dx);
            for (int i = 0; i < n; ++i){
                s.H(i, i) = hii;
            }
        }
        matrix_type tmp(n, n);
        linalg::gemm(1.0, s.H, linalg::trans(Z), 0.0, tmp);
        linalg::gemm(1.0, Z, tmp, 0.0, s.H);
        s.H += ss;
        s.dfx0 = s.dfx;
        return s;
    }

    ook::line_search::more_thuente search;
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
