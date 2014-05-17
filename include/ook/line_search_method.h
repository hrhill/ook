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

#ifndef OOK_LINE_SEARCH_METHOD_H_
#define OOK_LINE_SEARCH_METHOD_H_

#include <iostream>
#include <iomanip>
#include <boost/numeric/ublas/matrix.hpp>

#include "ook/norms.h"
#include "ook/message.h"
#include "ook/call_selector.h"

namespace ook{

template <typename Scheme>
struct line_search_method{

    typedef typename Scheme::state_type state_type;

    template <typename F, typename X, typename Options, typename Observer>
    std::tuple<ook::message, X>
    operator()(F obj_fun, X x, const Options& opts, Observer& observer)
    {
        typedef typename X::value_type real_type;
        typedef decltype(obj_fun(x)) result_type;
        typedef detail::call_selector<F, X, state_type,
                    std::tuple_size<result_type>::value> caller_type;

        const real_type epsilon = std::numeric_limits<real_type>::epsilon();
        const real_type dx_eps = sqrt(epsilon);
        const real_type df_eps = exp(log(epsilon)/real_type(3.0));

        state_type s = scheme.initialise(obj_fun, x);

        observer(s);
        while(true) {
            // Get descent direction and set up line search procedure.
            s = scheme.descent_direction(s);

            // Create line search function
            auto phi = [&s, &x, obj_fun](const real_type& a){
                caller_type::call(obj_fun, x + a * s.p, s);
                return std::make_pair(s.fx, inner_product(s.dfx, s.p));
            };

            // Store current fx value since line search overwrites the state values.
            const real_type fxk = s.fx;
            real_type dfx_dot_p = inner_product(s.dfx, s.p);
            std::tie(s.msg, s.a, s.fx, dfx_dot_p) = scheme.search(phi, s.fx, dfx_dot_p, 1.0, opts);

            if (s.msg != ook::message::convergence){
                break;
            }
            // check for warnings here too, perhaps contninuing if a descent direction can be found
            s.dx = s.a * s.p;
            x += s.dx;

            // Convergence criteria assessment base on p306 in Gill, Murray and Wright.
            const real_type theta = epsilon * (1.0 + fabs(s.fx));
            const bool u1 = (fxk - s.fx) <= theta;
            const bool u2 = ook::norm_infinity(s.dx) <=  dx_eps * (1.0 + ook::norm_infinity(x));
            const bool u3 = ook::norm_infinity(s.dfx) <= df_eps * (1.0 + fabs(s.fx));

            s.tag = detail::state_tag::iterate;
            observer(s);
            if ((u1 and u2) or u3){
                s.msg = ook::message::convergence;
                break;
            }
            s = scheme.update(s);
        }

        s.tag = detail::state_tag::final;
        observer(s);
        return std::make_pair(s.msg, x);
    }

    Scheme scheme;
};

} // ns ook

#endif
