
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

#include "ook/line_search/mcsrch.hpp"
#include "ook/line_search_method.hpp"

namespace ook{
/// \brief Implementation of the required steps of line_search_method
/// for Fletcher-Reeves method. The intitial step is a steepest descent
/// step. However, subsequent steps use the previous descent direction
/// in a linear combination with the steepest descent direction.
namespace beta{
struct fr
{
    double operator()(const vector& dxp, const vector& dxq, const vector& p) const
    {
        return (dxp, dxp) / (dxq, dxq);
    }
};

struct pr
{
    double operator()(const vector& dxp, const vector& dxq, const vector& p) const
    {
        return (dxp, dxp - dxq) / (dxq, dxq);
    }
};

struct hs
{
    double operator()(const vector& dxp, const vector& dxq, const vector& p) const
    {
        vector d = dxp - dxq;
        return (dxp, d) / (p, d);
    }
};

struct dy
{
    double operator()(const vector& dxp, const vector& dxq, const vector& p) const
    {
        return (dxp, dxp) / (p, dxp - dxq);
    }
};
} // ns beta

template <typename Beta>
struct nonlinear_cg_impl
{
    struct state
    {
        typedef vector vector_type;
        typedef double value_type;
        typedef matrix matrix_type;

    };

    /// \brief Initialize the scheme with \f$ \nabla f (x_0)\f$, and
    /// \f$ \beta = 0 \f$.
    template <typename State>
    explicit nonlinear_cg_impl(const State& s)
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
    vector
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
        beta = beta_fun(s.dfx, dfx, p);
        dfx = s.dfx;
    }

private:
    uint iter;
    vector dfx;
    vector p;
    double beta;
    Beta beta_fun;
};

/// \brief The Fletcher-Reeves algorithm.
template <typename F, typename Options, typename Observer>
typename line_search_method<nonlinear_cg_impl<beta::fr>, line_search::mcsrch>::state_type
fletcher_reeves(F f, const vector& x0, const Options& opts, Observer& observer)
{
    typedef nonlinear_cg_impl<beta::fr> scheme;
    line_search_method<scheme, line_search::mcsrch> method;
    return method(f, x0, opts, observer);
}

/// \brief The Polak-Ribiere algorithm.
template <typename F, typename Options, typename Observer>
typename line_search_method<nonlinear_cg_impl<beta::pr>, line_search::mcsrch>::state_type
polak_ribiere(F f, const vector& x0, const Options& opts, Observer& observer)
{
    typedef nonlinear_cg_impl<beta::pr> scheme;
    line_search_method<scheme, line_search::mcsrch> method;
    return method(f, x0, opts, observer);
}

/// \brief The Hestenes-Steifel algorithm.
template <typename F, typename Options, typename Observer>
typename line_search_method<nonlinear_cg_impl<beta::hs>, line_search::mcsrch>::state_type
hestenes_steifel(F f, const vector& x0, const Options& opts, Observer& observer)
{
    typedef nonlinear_cg_impl<beta::hs> scheme;
    line_search_method<scheme, line_search::mcsrch> method;
    return method(f, x0, opts, observer);
}

/// \brief The Dai-Yuan algorithm.
template <typename F, typename Options, typename Observer>
typename line_search_method<nonlinear_cg_impl<beta::dy>, line_search::mcsrch>::state_type
dai_yuan(F f, const vector& x0, const Options& opts, Observer& observer)
{
    typedef nonlinear_cg_impl<beta::dy> scheme;
    line_search_method<scheme, line_search::mcsrch> method;
    return method(f, x0, opts, observer);
}

} //ns ook

#endif
