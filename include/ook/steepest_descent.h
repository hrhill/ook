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

#ifndef OOK_LINE_SEARCH_METHODS_STEEPEST_DESCENT_H_
#define OOK_LINE_SEARCH_METHODS_STEEPEST_DESCENT_H_

#include <tuple>

#include "ook/state.h"
#include "ook/line_search_method.h"
#include "ook/line_search/more_thuente.h"

namespace ook{
namespace detail{

/// \brief Implementation of the required steps of line_search_method
/// for the steepes descent method.
template <typename X>
struct steepest_descent
{
    typedef X vector_type;
    typedef typename X::value_type value_type;
    typedef state<X> state_type;

    template <typename F>
    state_type
    initialise(F objective_function, const X& x0)
    {
        state_type s(x0.size());
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
    update(state_type s){
        return s;
    }
};

} // ns detail

/// \brief The Steepest descent algorithm.
template <typename F, typename X, typename Options, typename Observer>
std::tuple<ook::message, X>
steepest_descent(F obj_fun, const X& x0, const Options& opts, Observer& observer)
{
    typedef detail::steepest_descent<X> scheme;
    line_search_method<scheme, ook::line_search::more_thuente> method;
    return method(obj_fun, x0, opts, observer);
}

} //ns ook


#endif
