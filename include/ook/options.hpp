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

#ifndef OOK_OPTIONS_HPP_
#define OOK_OPTIONS_HPP_

#include <cassert>
#include <limits>

namespace ook{

template <typename T>
struct options
{
    options()
    :
        ftol(T(1e-03)),
        gtol(T(9e-01)),
        stpmin(std::numeric_limits<T>::epsilon()),
        stpmax(T(5.0)),
        max_function_evaluations(2000)
    {}

    options(T ftol, T gtol, T stpmin, T stpmax)
    :
        ftol(ftol), gtol(gtol),  stpmin(stpmin), stpmax(stpmax)
    {
        assert(ftol > T(0) && "ftol <= 0.0");
        assert(gtol > T(0) && "gtol <= 0.0");
        assert(stpmin >= T(0.0) && "stpmin < 0.0");
        assert(stpmax > stpmin && "stpmax <= stpmin.");
    }

    // tolerance for the sufficient decrease condition.
    T ftol;
    // tolerance for the curvature condition.
    T gtol;
    // lower bound for the step.
    T stpmin;
    // upper bound for the step.
    T stpmax;

    uint max_function_evaluations;
};

} // ns ook

#endif
