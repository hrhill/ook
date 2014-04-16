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
#include <boost/numeric/ublas/matrix.hpp>

#include "ook/norms.h"
#include "ook/state.h"
#include "ook/line_search_method.h"
#include "ook/line_search/more_thuente.h"

namespace ook{
namespace detail{

/// \brief Implementation of the required steps of line_search_method
/// for BFGS method.
template <typename X>
struct bfgs{
    typedef X vector_type;
    typedef typename X::value_type value_type;
    typedef state<X> state_type;

    template <typename F>
    static
    state_type
    initialise(F objective_function, const X& x0)
    {
        state_type s(x0.size(), true);
        std::tie(s.fx, s.dfx) = objective_function(x0);
        s.dfx0 = s.dfx;
        return s;
    }

    static
    vector_type
    descent_direction(state_type& s)
    {
        ++s.iteration;
        s.p = -boost::numeric::ublas::prod(s.H, s.dfx);
        return s.p;
    }

    static
    state_type
    update(state_type s){
        namespace ublas = boost::numeric::ublas;
        typedef ublas::matrix<double, ublas::column_major> matrix_type;

        X dx = s.a * s.p;
        X y(s.dfx - s.dfx0);
        const int n = s.dfx.size();

        const value_type rho = 1.0/inner_product(y, dx);
        matrix_type Z(ublas::identity_matrix<double>(n) - rho * ublas::outer_prod(dx, y));
        matrix_type ss = rho * ublas::outer_prod(dx, dx);

        if (s.iteration == 1){
            const value_type hii = inner_product(dx, dx);
            for (int i = 0; i < n; ++i){
                s.H(i, i) = hii;
            }
        }
        s.H = ublas::prod(Z, matrix_type(ublas::prod(s.H, ublas::trans(Z)))) + ss;
        s.dfx0 = s.dfx;
        return s;
    }
};

} // ns detail

/// \brief The Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm.
template <typename F, typename X, typename Options, typename Observer>
std::tuple<ook::message, X>
bfgs(F obj_fun, const X& x0, const Options& opts, Observer& observer)
{
    typedef detail::bfgs<X> scheme;
    line_search_method<scheme, ook::line_search::more_thuente> method;
    return method.run(obj_fun, x0, opts, observer);
}

} //ns ook


#endif
