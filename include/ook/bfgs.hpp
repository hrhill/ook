
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

#ifndef OOK_LINE_SEARCH_METHODS_BFGS_HPP_
#define OOK_LINE_SEARCH_METHODS_BFGS_HPP_

#include <tuple>

#include "ook/matrix.hpp"
#include "ook/vector.hpp"

#include "ook/line_search/mcsrch.hpp"
#include "ook/line_search_method.hpp"

namespace ook
{
/// \brief Implementation of the required steps of line_search_method
/// for BFGS method.
struct bfgs_impl
{
    struct state
    {
        matrix H;
    };

    template <typename State>
    explicit bfgs_impl(const State& s)
        : dfx(s.dfx), p(s.dfx.size()), H(s.dfx.size(), s.dfx.size(), 0.0)
    {
        const int n = dfx.size();
        for (int i = 0; i < n; ++i)
        {
            H(i, i) = 1.0;
        }
    }

    template <typename State>
    vector
    descent_direction(const State& s)
    {
        p = -H * s.dfx;
        return p;
    }

    template <typename State>
    void
    update(const State& s)
    {
        vector yk = s.dfx - dfx;
        const double rho = 1.0 / (yk, s.dx);
        const int n = s.dfx.size();
        matrix Z = -rho * (s.dx * trans(yk));

        for (int i = 0; i < n; ++i)
            Z(i, i) += 1.0;

        matrix tmp = H * trans(Z);
        H = Z * tmp;
        matrix ss = rho * s.dx * trans(s.dx);
        H += ss;
        dfx = s.dfx;
    }

private:
    vector dfx;
    vector p;
    matrix H;
};

/// \brief The Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm.
template <typename F, typename Options, typename Observer>
typename line_search_method<bfgs_impl, line_search::mcsrch>::state_type
bfgs(F f, const vector& x0, const Options& opts, Observer& observer)
{
    line_search_method<bfgs_impl, line_search::mcsrch> method;
    return method(f, x0, opts, observer);
}

} // ns ook

#endif
