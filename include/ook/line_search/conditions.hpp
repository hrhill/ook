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

#ifndef OOK_LINE_SEARCH_CONDITIONS_HPP_
#define OOK_LINE_SEARCH_CONDITIONS_HPP_

#include <cassert>
#include <cmath>

namespace ook{

namespace line_search{

template <typename T>
bool
sufficient_decrease_condition(const T& fxap,
                              const T& fx,
                              const T& c,
                              const T& a,
                              const T& dfx_dot_p)
{
    assert(T(0.0) < c && c < T(1.0));
    return fxap <= fx + c * a * dfx_dot_p;
}

template <typename T>
bool
curvature_condition(const T& dfxap_dot_p,
                    const T& c,
                    const T& dfx_dot_p)
{
    assert(T(0.0) < c && c < T(1.0));
    return fabs(dfxap_dot_p) <= c * fabs(dfx_dot_p);
}

template <typename T>
bool
strong_wolfe_conditions(const T& fx,
                        const T& fxap,
                        const T& dfx_dot_p,
                        const T& dfxap_dot_p,
                        const T& a,
                        const T& c1 = T(1e-04),
                        const T& c2 = T(9e-01))
{
    return sufficient_decrease_condition(fxap, fx, c1, a, dfx_dot_p)
            && curvature_condition(dfx_dot_p, dfxap_dot_p, c2);
}

} // ns line_search

} // ns ook

#endif
