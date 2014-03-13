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

#ifndef OOK_LINE_SEARCH_MORE_THUENTE_H_
#define OOK_LINE_SEARCH_MORE_THUENTE_H_

#include <tuple>
#include <exception>
#include <cassert>

#include "ook/message.h"
#include "ook/line_search/conditions.h"
#include "ook/line_search/more_thuente/search.h"

namespace ook{
namespace line_search{

struct more_thuente{

    template <typename F, typename T, typename Options>
    static
    std::tuple<message, T, T, T>
    search(F phi, T phi0, T dphi0, T a, const Options& opts)
    {
        if(a <= opts.stpmin)
            return std::make_tuple(
                        message::warning_stp_eq_stpmin, 0, phi0, dphi0);
        if(a >= opts.stpmax)
            return std::make_tuple(
                        message::warning_stp_eq_stpmax, 0, phi0, dphi0);
        if (dphi0 > T(0.0)) {
            return std::make_tuple(
                        message::search_direction_is_not_a_descent_direction,
                         0, phi0, dphi0);
        }

        T phia, dphia;
        std::tie(phia, dphia) = phi(a);
        const T width = opts.stpmax - opts.stpmin;
        detail::search<T> search(0.0, phi0, dphi0, a, phia, dphia, 0.0, 5.0 * a, width, 2.0 * width);
        message msg;
        while(true){
            std::tie(msg, a) = search(a, phia, dphia, opts);
            if (msg != message::update) break;
            std::tie(phia, dphia) = phi(a);
        }
        return std::make_tuple(msg, a, phia, dphia);
    }
};

} // ns linesearch
} // ns ook

#endif
