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

#ifndef OOK_LINE_SEARCH_METHOD_HPP_
#define OOK_LINE_SEARCH_METHOD_HPP_

#include <iostream>
#include <iomanip>

#include "linalg.hpp"

#include "ook/message.hpp"
#include "ook/state.hpp"
#include "ook/call_selector.hpp"

namespace ook{

template <typename X>
struct lsm_state
{
    typedef typename linalg::associated_matrix<X>::type matrix_type;

    lsm_state(const X& x)
    :
        iteration(0),
        nfev(0),
        a(0),
        fx(0),
        gnorm(0),
        xnorm(0),
        tag(state_tag::init),
        x(x),
        dfx(x.size(), 0),
        dx(x.size(), 0),
        H(x.size(), x.size())
    {

    }

    uint iteration;
    uint nfev;
    double a;
    double fx;
    double gnorm;
    double xnorm;
    state_tag tag;
    message msg;
    X x;
    X dfx;
    X dx;
    matrix_type H;
};

template <typename Scheme>
struct line_search_method
{

    template <typename F, typename X, typename Options, typename Observer>
    std::tuple<ook::message, X>
    operator()(F obj_fun, X x, const Options& opts, Observer& observer)
    {
        typedef detail::call_selector<
                    std::tuple_size<decltype(obj_fun(x))>::value> caller_type;
        typedef typename X::value_type real_type;

        const real_type epsilon = std::numeric_limits<real_type>::epsilon();
        const real_type dx_eps = sqrt(epsilon);
        const real_type df_eps = exp(log(epsilon)/real_type(3.0));

        lsm_state<X> state(x);
        state = caller_type::call(obj_fun, x, state);
        Scheme scheme(state);

        observer(state);

        while(true)
        {
            // Get descent direction and set up line search procedure.
            X p = scheme.descent_direction(state);

            // Create line search function
            auto phi = [&p, &x, obj_fun, &state](real_type a)
            {
                ++state.nfev;
                state = caller_type::call(obj_fun, x + a * p, state);
                return std::make_pair(state.fx, linalg::inner_prod(state.dfx, p));
            };

            // Store current fx value since line search overwrites the state values.
            const real_type fxk = state.fx;
            real_type dfx_dot_p = linalg::inner_prod(state.dfx, p);
            std::tie(state.msg, state.a, state.fx, dfx_dot_p)
                = scheme.search(phi, state.fx, dfx_dot_p, 1.0, opts);

            if (state.msg != ook::message::convergence){
                break;
            }

            state.dx = state.a * p;
            x += state.dx;
            state.xnorm = linalg::norm_infinity(state.dx);
            state.gnorm = linalg::norm_infinity(state.dfx);

            // Convergence criteria assessment base on p306 in Gill, Murray and Wright.
            const real_type theta = epsilon * (1.0 + fabs(state.fx));
            const bool u1 = (fxk - state.fx) <= theta;
            const bool u2 = state.xnorm <=  dx_eps * (1.0 + linalg::norm_infinity(x));
            const bool u3 = state.gnorm <= df_eps * (1.0 + fabs(state.fx));

            state.tag = state_tag::iterate;
            observer(state);
            if ((u1 and u2) or u3){
                state.msg = ook::message::convergence;
                break;
            }
            scheme.update(state);
            ++state.iteration;
        }
        state.tag = state_tag::final;
        observer(state);
        return std::make_pair(state.msg, x);
    }
};

} // ns ook

#endif
