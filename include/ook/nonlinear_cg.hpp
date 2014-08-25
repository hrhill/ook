
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
namespace beta{
struct fr
{
    template <typename X>
    auto operator()(const X& dxp, const X& dxq, const X& p) const
    -> remove_const_reference<decltype(dxp[0])>
    {
        return linalg::inner_prod(dxp, dxp)/linalg::inner_prod(dxq, dxq);
    }
};

struct pr
{
    template <typename X>
    auto operator()(const X& dxp, const X& dxq, const X& p) const
    -> remove_const_reference<decltype(dxp[0])>
    {
        return linalg::inner_prod(dxp,
                static_cast<const X&>(dxp - dxq))/linalg::inner_prod(dxq, dxq);
    }
};

struct hs
{
    template <typename X>
    auto operator()(const X& dxp, const X& dxq, const X& p) const
    -> remove_const_reference<decltype(dxp[0])>
    {
        X d = dxp - dxq;
        return linalg::inner_prod(dxp, d)/linalg::inner_prod(p, d);
    }
};

struct dy
{
    template <typename X>
    auto operator()(const X& dxp, const X& dxq, const X& p) const
    -> remove_const_reference<decltype(dxp[0])>
    {
        return linalg::inner_prod(dxp, dxp)/linalg::inner_prod(p,
                    static_cast<const X&>(dxp - dxq));
    }
};
} // ns beta

template <typename X, typename Beta>
struct nonlinear_cg_impl
{
    typedef X vector_type;
    typedef typename std::remove_reference<decltype(X()[0])>::type value_type;

    struct state
    {
        typedef X vector_type;
        typedef typename std::remove_reference<decltype(X()[0])>::type value_type;

        value_type fx;
        vector_type dfx;
    };

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
    /// \f[ p_{k+1} = -\nabla f(x_k) + \beta_k p_k\f]
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

    /// \brief Update \f$\beta\f$ and \f$\nabla f(x)\f$.
    /// \tparam State The obtimiser state type.
    /// \param s The current state.
    template <typename State>
    void
    update(const State& s)
    {
        const int n = dfx.size();
        beta = beta_fun(s.dfx, dfx, p);
        dfx = s.dfx;
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
fletcher_reeves(F f, const X& x0, const Options& opts, Observer& observer)
{
    typedef nonlinear_cg_impl<X, beta::fr> scheme;
    line_search_method<scheme> method;
    return method(f, x0, opts, observer);
}

/// \brief The Polak-Ribiere algorithm.
template <typename F, typename X, typename Options, typename Observer>
std::tuple<ook::message, X>
polak_ribiere(F f, const X& x0, const Options& opts, Observer& observer)
{
    typedef nonlinear_cg_impl<X, beta::pr> scheme;
    line_search_method<scheme> method;
    return method(f, x0, opts, observer);
}

/// \brief The Hestenes-Steifel algorithm.
template <typename F, typename X, typename Options, typename Observer>
std::tuple<ook::message, X>
hestenes_steifel(F f, const X& x0, const Options& opts, Observer& observer)
{
    typedef nonlinear_cg_impl<X, beta::hs> scheme;
    line_search_method<scheme> method;
    return method(f, x0, opts, observer);
}

/// \brief The Dai-Yuan algorithm.
template <typename F, typename X, typename Options, typename Observer>
std::tuple<ook::message, X>
dai_yuan(F f, const X& x0, const Options& opts, Observer& observer)
{
    typedef nonlinear_cg_impl<X, beta::dy> scheme;
    line_search_method<scheme> method;
    return method(f, x0, opts, observer);
}

} //ns ook

#endif
