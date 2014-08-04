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

#ifndef OOK_LINE_SEARCH_METHODS_MESSAGE_HPP_
#define OOK_LINE_SEARCH_METHODS_MESSAGE_HPP_

#include <iostream>
#include <type_traits>

namespace ook{

enum class message{
    start,
    update,
    warning_rounding_error_prevents_progress,
    warning_xtol_satisfied,
    warning_stp_eq_stpmax,
    warning_stp_eq_stpmin,
    warning_max_line_search_attempts_reached,
    convergence,
    search_direction_is_not_a_descent_direction,
    error_step_less_than_stpmin,
    error_step_greater_than_stpmax
};

const char* message_string[] = {"start",
    "update",
    "warning_rounding_error_prevents_progress",
    "warning_xtol_satistfied",
    "warning_stp_eq_stpmax",
    "warning_stp_eq_stpmin",
    "warning_max_line_search_attempts_reached",
    "convergence",
    "search_direction_is_not_a_descent_direction",
    "error_step_less_than_stpmin",
    "error_step_greater_than_stpmax"
};

std::ostream&
operator<<(std::ostream& os, const message& tv)
{
    auto id = static_cast<std::underlying_type<message>::type>(tv);
    return os << std::string(message_string[id]);
}

} // ns ook

#endif
