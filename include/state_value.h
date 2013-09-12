#ifndef OOK_STATE_VALUES_H_
#define OOK_STATE_VALUES_H_

namespace ook{

enum class state_value{
    start,
    update,
    warning_rounding_error_prevents_progress,
    warning_xtol_satisfied,
    warning_stp_eq_stpmax,
    warning_stp_eq_stpmin,
    convergence
};

const char* state_value_string[] = {"start",
    "update",
    "warning_rounding_error_prevents_progress",
    "warning_xtol_satistfied",
    "warning_stp_eq_stpmax",
    "warning_stp_eq_stpmin",
    "convergence"}; 

std::ostream&
operator<<(std::ostream& os, const state_value& tv)
{
    auto id = static_cast<std::underlying_type<state_value>::type>(tv);
    return os << std::string(state_value_string[id]);
}

} // ns ook

#endif