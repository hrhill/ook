
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
#include <type_traits>

#include "ook/vector.hpp"

#include "ook/message.hpp"
#include "ook/call_selector.hpp"
#include "ook/line_search/mcsrch.hpp"

namespace ook{

template <typename SchemeState>
struct lsm_state : public SchemeState
{
    /// \brief State for use with line search method.
    enum class tag{init, iterate, final};

    lsm_state()
    :
        iteration(0),
        nfev(0),
        a(0),
        fx(0),
        gnorm(0),
        xnorm(0),
        state(lsm_state::tag::init)
    {}

    friend
    std::ostream&
    operator<<(std::ostream& out, const lsm_state& s)
    {
        if (s.state == lsm_state::tag::init)
        {
            out << "\n"
                << std::setw(8) << "n"
                << std::setw(8) << "nfev"
                << std::scientific
                << std::setw(16) << "a"
                << std::setw(16) << "fx"
                << std::setw(16) << "max ||dfx||"
                << std::setw(16) << "max ||dx||";
        }

        if (s.state == lsm_state::tag::iterate)
        {
            out << std::setw(8) << s.iteration
                << std::setw(8) << s.nfev
                << std::scientific
                << std::setw(16) << s.a
                << std::setw(16) << s.fx
                << std::setw(16) << s.gnorm
                << std::setw(16) << s.xnorm;
        }

        if (s.state == lsm_state::tag::final)
        {
            out << "\nstatus : " << s.msg
                << "\n"
                << std::setw(8) << "iter"
                << std::setw(8) << "nfev"
                << std::setw(16) << "fx"
                << std::setw(16) << "max ||dfx||"
                << std::setw(16) << "max ||dx||\n"
                << std::setw(8) << s.iteration
                << std::setw(8) << s.nfev
                << std::scientific
                << std::setw(16) << s.fx
                << std::setw(16) << s.gnorm
                << std::setw(16) << s.xnorm;
        }
        return out;
    }

    uint iteration;
    uint nfev;
    double a;
    double fx;
    double gnorm;
    double xnorm;
    tag state;
    message msg;
    vector x;
    vector dfx;
    vector dx;
};

template <typename Scheme, typename LineSearch>
struct line_search_method
{
    typedef lsm_state<typename Scheme::state> state_type;

    template <typename F, typename X, typename Options, typename Observer>
    state_type
    operator()(F obj_fun, X x, const Options& opts, Observer& observer) const
    {
        LineSearch line_search;

        typedef detail::call_selector<
                std::tuple_size<decltype(obj_fun(x))>::value> caller_type;

        const double epsilon = std::numeric_limits<double>::epsilon();
        const double dx_eps = sqrt(epsilon);
        const double df_eps = exp(log(epsilon)/double(3.0));

        state_type state = caller_type::call(obj_fun, x, state_type());
        Scheme scheme(state);
        observer(state);

        while(true)
        {
            // Get descent direction and set up line search procedure.
            X p = scheme.descent_direction(state);

            // Create line search function
            auto phi = [&p, &x, obj_fun, &state](double a)
            {
                ++state.nfev;
                state = caller_type::call(obj_fun, x + a * p, state);
                return std::make_pair(state.fx, (state.dfx, p));
            };

            // Store current fx value since line search overwrites the state
            // values.
            const double fxk = state.fx;
            double dfx_dot_p = (state.dfx, p);

            if (dfx_dot_p >= 0.0){
                state.msg =
                    message::search_direction_is_not_a_descent_direction;
                break;
            }

            std::tie(state.msg, state.a, state.fx, dfx_dot_p)
                = line_search(phi, state.fx, dfx_dot_p, 1.0, opts);

            if (state.msg != message::convergence) break;

            state.dx = state.a * p;
            x += state.dx;
            state.xnorm = norm_inf(state.dx);
            state.gnorm = norm_inf(state.dfx);

            // Convergence criteria assessment base on p306 in Gill, Murray
            // and Wright.
            const double theta = epsilon * (1.0 + fabs(state.fx));
            const bool u1 = (fxk - state.fx) <= theta;
            const bool u2 = state.xnorm
                            <= dx_eps * (1.0 + norm_inf(x));
            const bool u3 = state.gnorm <= df_eps * (1.0 + fabs(state.fx));

            state.state = state_type::tag::iterate;
            observer(state);
            if ((u1 and u2) or u3){
                state.msg = message::convergence;
                break;
            }
            if (state.iteration >= opts.maxiter){
                state.msg = message::maximum_iterations_reached;
                break;
            }

            scheme.update(state);
            ++state.iteration;
        }
        state.state = state_type::tag::final;
        state.x = x;
        observer(state);
        return state;
    }
};

} // ns ook

#endif
