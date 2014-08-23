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

#ifndef OOK_LINE_SEARCH_METHODS_NONLINEAR_CG_HPP_
#define OOK_LINE_SEARCH_METHODS_NONLINEAR_CG_HPP_

#include <tuple>
#include <algorithm>
#include <type_traits>

#include "linalg.hpp"

#include "ook/line_search_method.hpp"

namespace ook{
/// \brief Implementation of the required steps of line_search_method
/// for Fletcher-Reeves method. The intitial step is a steepest descent
/// step. However, subsequent steps use the previous descent direction
/// in a linear combination with the steepest descent direction.

struct beta_fr
{
    template <typename X>
    auto operator()(const X& dxp, const X& dxq, const X& p)
    {
        return linalg::inner_prod(dxp, dxp)/linalg::inner_prod(dxq, dxq);
    }
};

struct beta_pr
{
    template <typename X>
    auto operator()(const X& dxp, const X& dxq, const X& p)
    {
        return linalg::inner_prod(dxp, dxp - dxq)/linalg::inner_prod(dxq, dxq);
    }
};

struct beta_hs
{
    template <typename X>
    auto operator()(const X& dxp, const X& dxq, const X& p)
    {
        return linalg::inner_prod(dxp, dxp - dxq)/linalg::inner_prod(p, dxp - dxq);
    }
};

struct beta_dy
{
    template <typename X>
    auto operator()(const X& dxp, const X& dxq, const X& p)
    {
        return linalg::inner_prod(dxp, dxp)/linalg::inner_prod(p, dxp - dxq);
    }
};


template <typename X, typename Beta>
struct nonlinear_cg_impl
{
    typedef X vector_type;
    typedef typename std::remove_reference<decltype(X()[0])>::type value_type;

    /// \brief Initialize the scheme with \f$ \nabla f (x_0)\f$, and
    /// \f$ \beta = 0 \f$.
    template <typename State>
    nonlinear_cg_impl(const State& s)
    :
        iter(0),
        dfx(s.dfx),
        p(s.dfx.size(), 0),
        beta(0)
    {}

    /// \brief Calculate and return the descent direction given by
    ///
    /// \tparam State The obtimiser state type.
    /// \param s The current state.
    /// \return The descent direction.
    template <typename State>
    X
    descent_direction(const State& s)
    {
        p = -s.dfx + beta * p;
        return p;
    }

    /// \brief Update \f$\beta\f$, reseting after n steps.
    template <typename State>
    void
    update(const State& s)
    {
        const int n = dfx.size();
        beta = (iter % n) ? beta_fun(s.dfx, dfx, p) : 0;
        dfx = s.dfx;
        ++iter;
    }

private:
    uint iter;
    vector_type dfx;
    vector_type p;
    value_type beta;
    Beta beta_fun;
};

/// \brief The Fletcher-Reeves algorithm.
template <typename F, typename X, typename Options, typename Observer>
std::tuple<ook::message, X>
nonlinear_cg(F f, const X& x0, const Options& opts, Observer& observer)
{
    typedef nonlinear_cg_impl<X> scheme;
    line_search_method<scheme> method;
    return method(f, x0, opts, observer);
}
} //ns ook

#endif
