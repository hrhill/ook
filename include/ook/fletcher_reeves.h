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

#ifndef OOK_LINE_SEARCH_METHODS_FLETCHER_REEVES_H_
#define OOK_LINE_SEARCH_METHODS_FLETCHER_REEVES_H_

#include <tuple>
#include <algorithm>

#include "ook/state.h"
#include "ook/line_search_method.h"
#include "ook/line_search/more_thuente.h"

namespace ook{
namespace detail{

template <typename X>
struct
fr_state{
    typedef X vector_type;
    typedef typename X::value_type value_type;

    fr_state(const int n = 0)
    :
        fx(0),
        dfx(n),
        dfx0(n),
        p(n),
        dx(n),
        a(1),
        beta(0),
        iteration(0),
        nfev(0),
        tag(state_tag::init)
    {}

    value_type fx;
    vector_type dfx;
    vector_type dfx0;
    vector_type p;
    vector_type dx;
    value_type a;
    value_type beta;
    int iteration;
    int nfev;
    state_tag tag;
    message msg;
};

/// \brief Implementation of the required steps of line_search_method
/// for Fletcher-Reeves method.
template <typename X>
struct fletcher_reeves{
    typedef X vector_type;
    typedef typename X::value_type value_type;
    typedef fr_state<X> state_type;

    template <typename F>
    state_type
    initialise(F objective_function, const X& x0)
    {
        state_type s(x0.size());
        std::tie(s.fx, s.dfx) = objective_function(x0);
        s.dfx0 = s.dfx;
        return s;
    }

    state_type
    descent_direction(state_type s)
    {
        ++s.iteration;
        s.p = -s.dfx + s.beta * s.p;
        return s;
    }

    state_type
    update(state_type s)
    {
        s.beta = inner_product(s.dfx, s.dfx)/inner_product(s.dfx0, s.dfx0);
        s.dfx0 = s.dfx;
        return s;
    }

    ook::line_search::more_thuente search;
};

} // ns detail

/// \brief The Fletcher-Reeves algorithm.
template <typename F, typename X, typename Options, typename Observer>
std::tuple<ook::message, X>
fletcher_reeves(F f, const X& x0, const Options& opts, Observer& observer)
{
    typedef detail::fletcher_reeves<X> scheme;
    line_search_method<scheme> method;
    return method(f, x0, opts, observer);
}

} //ns ook

#endif
