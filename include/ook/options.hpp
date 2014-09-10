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

#include "ook/line_search/options.hpp"

namespace ook{

template <typename T>
struct options : line_search::options<T>
{
    options(T eps = T(1e-05), uint maxiter = 2000)
    :
        eps(eps), maxiter(maxiter)
    {}

    // tolerance for the first order convergence criteria
    T eps;
    // maximum number of iterations
    uint maxiter;
};

} // ns ook

#endif
