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

#ifndef OOK_LINE_SEARCH_METHODS_STEEPEST_DESCENT_HPP_
#define OOK_LINE_SEARCH_METHODS_STEEPEST_DESCENT_HPP_

#include <tuple>
#include <type_traits>

#include "linalg/operations.hpp"

#include "ook/state.hpp"
#include "ook/line_search_method.hpp"
#include "ook/line_search/more_thuente.hpp"

namespace ook{

namespace detail{

template <typename X>
struct
sd_state{
    typedef X vector_type;
    typedef typename std::remove_reference<decltype(X()[0])>::type value_type;

    sd_state(const int n = 0)
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
/// for the steepes descent method.
template <typename X>
struct steepest_descent
{
    typedef X vector_type;
    typedef typename std::remove_reference<decltype(X()[0])>::type value_type;
    typedef sd_state<X> state_type;

    template <typename F>
    state_type
    initialise(F objective_function, const X& x0)
    {
        state_type s(linalg::size(x0));
        std::tie(s.fx, s.dfx) = objective_function(x0);
        return s;
    }

    state_type
    descent_direction(state_type s)
    {
        ++s.iteration;
        s.p = -s.dfx;
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

/// \brief The Steepest descent algorithm.
template <typename F, typename X, typename Options, typename Observer>
std::tuple<ook::message, X>
steepest_descent(F f, const X& x0, const Options& opts, Observer& observer)
{
    typedef detail::steepest_descent<X> scheme;
    line_search_method<scheme> method;
    return method(f, x0, opts, observer);
}

} //ns ook

#endif
