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

#ifndef OOK_LINE_SEARCH_METHODS_FLETCHER_REEVES_HPP_
#define OOK_LINE_SEARCH_METHODS_FLETCHER_REEVES_HPP_

#include <tuple>
#include <algorithm>
#include <type_traits>

#include "linalg.hpp"

#include "ook/line_search_method.hpp"
#include "ook/line_search/more_thuente.hpp"

namespace ook{
namespace detail{

template <typename X>
struct
fletcher_reeves_state
{
    typedef X vector_type;
    typedef typename std::remove_reference<decltype(X()[0])>::type value_type;

    fletcher_reeves_state(const X& dfx)
    :
        dfx(dfx),
        p(dfx.size(), 0),
        beta(0)
    {}

    vector_type dfx;
    vector_type p;
    value_type beta;
};

/// \brief Implementation of the required steps of line_search_method
/// for Fletcher-Reeves method. The intitial step is a steepest descent
/// step. However, subsequent steps use the previous descent direction
/// is a linear combination with the steepest descent direction.
template <typename X>
struct fletcher_reeves
{
    typedef fletcher_reeves_state<X> state_type;

    template <typename State>
    fletcher_reeves(const State& s)
    :
        state(s.dfx)
    {}

    template <typename State>
    X
    descent_direction(const State& s)
    {
        state.p = -s.dfx + state.beta * state.p;
        return state.p;
    }

    template <typename State>
    void
    update(const State& s)
    {
        state.beta = linalg::inner_prod(s.dfx, s.dfx)/linalg::inner_prod(state.dfx, state.dfx);
        state.dfx = s.dfx;
    }

    state_type state;
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
