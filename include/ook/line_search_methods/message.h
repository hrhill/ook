#ifndef OOK_LINE_SEARCH_METHODS_MESSAGE_H_
#define OOK_LINE_SEARCH_METHODS_MESSAGE_H_

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
    search_direction_is_not_a_descent_direction
};

const char* message_string[] = {"start",
    "update",
    "warning_rounding_error_prevents_progress",
    "warning_xtol_satistfied",
    "warning_stp_eq_stpmax",
    "warning_stp_eq_stpmin",
    "warning_max_line_search_attempts_reached",
    "convergence",
    "search_direction_is_not_a_descent_direction"}; 

std::ostream&
operator<<(std::ostream& os, const message& tv)
{
    auto id = static_cast<std::underlying_type<message>::type>(tv);
    return os << std::string(message_string[id]);
}

} // ns ook

#endif