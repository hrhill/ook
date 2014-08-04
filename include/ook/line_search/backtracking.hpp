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

#ifndef OOK_LINE_SEARCH_BACKTRACKING_HPP_
#define OOK_LINE_SEARCH_BACKTRACKING_HPP_

#include <tuple>
#include <exception>
#include <cassert>

#include "ook/message.hpp"
#include "ook/line_search/conditions.hpp"

namespace ook{
namespace line_search{

struct backtracking{
    template <typename F, typename T, typename Options>
    std::tuple<message, T, T, T>
    operator()(F phi, T phi0, T dphi0, T a, const Options& opts)
    {
        T phia, dphia;
        T rho(0.9); // Need to pass this as an option.
        ook::message msg;
        while(true){
            if(a < opts.stpmin){
                msg = message::warning_stp_eq_stpmin;
                break;
            }
            if(a > opts.stpmax){
                msg = message::warning_stp_eq_stpmax;
                break;
            }

            std::tie(phia, dphia) = phi(a);
            if (sufficient_decrease_condition(phia, phi0, opts.ftol, a, dphi0)){
                msg = message::convergence;
                break;
            }
            a *= rho;
        }
        return std::make_tuple(msg, a, phia, dphia);
    }
};


} // ns line_search
} // ns ook

#endif
