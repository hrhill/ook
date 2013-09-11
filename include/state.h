#ifndef OOK_STATE_H_
#define OOK_STATE_H_

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>

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

struct state{
    typedef double real_type;

    state(double fx_, double fxap_, double dfx_dot_p_, double dfxap_dot_p_, double a_, state_value value_, int nfev_ = 0)
    : 
        fx(fx_), 
        fxap(fxap_), 
        dfx_dot_p(dfx_dot_p_),
        dfxap_dot_p(dfxap_dot_p_),
        a(a_),
        value(value_),
        nfev(nfev_)
    {}

    double fx;
    double fxap;
    double dfx_dot_p;
    double dfxap_dot_p;
    double a;
    state_value value;
    int nfev;
    boost::numeric::ublas::vector<double> x;
    boost::numeric::ublas::vector<double> p;
};

}  // ns ook

#endif