#include <iostream>
#include <algorithm>
#include <tuple>
#include <iomanip>

#include <boost/math/constants/constants.hpp>

#include "line_search.h"

template <typename Enumeration>
auto as_integer(Enumeration const value)
    -> typename std::underlying_type<Enumeration>::type
{
    return static_cast<typename std::underlying_type<Enumeration>::type>(value);
}

const char* task_value_string[] = {"start",
    "fg",
    "warning_rounding_error_prevents_progress",
    "warning_xtol_satistfied",
    "warning_stp_eq_stpmax",
    "warning_stp_eq_stpmin",
    "CONVERGENCE"}; // for testing, don't want to change this yet

std::ostream&
operator <<(std::ostream& os, const task_value& tv)
{
    auto id = static_cast<std::underlying_type<task_value>::type>(tv);
   os << std::string(task_value_string[id]);
   return os;
}

std::tuple<double, double>
phi51(double a, double b){
    return std::make_tuple(-a/(a * a + b),
                          (a * a - b)/((a * a + b) * (a * a + b)));
}

std::tuple<double, double>
phi52(double a, double b){
    const double apb = (a + b);
    return std::make_pair(std::pow(apb, 5) - 2 * std::pow(apb, 4),
                          5.0 * std::pow(apb, 4) - 8 * std::pow(apb, 3));
}

std::tuple<double, double>
phi53(double a, double b, double c){
    const double pi = boost::math::constants::pi<double>();

    double phi0;
    double dphi0;
    if (a <= 1 - b){
        phi0 = 1 - a;
        dphi0 = -1;
    }else if (a >= 1 + b){
        phi0 = a - 1;
        dphi0 = 1;        
    }else{
        phi0 = 0.5 * std::pow(a - 1, 2)/b + 0.5 * b;
        dphi0 = (a - 1)/b;        
    }
    return std::make_tuple(phi0 + 2 * (1 - b)/(c * pi) * sin(0.5 * c * pi * a), 
                           dphi0 + (1 - b) * cos(0.5 * c * pi * a));
}

std::tuple<double, double>
phi54(double a, double b1, double b2){

    auto gamma = [](const double b) -> double { 
        return sqrt(1.0 + std::pow(b, 2)) - b;
    };

    const double t1 = sqrt(std::pow(1 - a, 2) + pow(b2, 2));
    const double t2 = sqrt(std::pow(a, 2) + std::pow(b1, 2));
    return std::make_tuple(
            gamma(b1) * t1 + gamma(b2) * t2,
          - gamma(b1) * (1 - a)/t1 + gamma(b2) * a/t2);
}

template <typename F>
void
more_thuente_line_search(F phi, const double stp0, const double& mu, const double eta){

    options opts{mu, 
                 eta,
                 std::numeric_limits<double>::epsilon(), 
                 0.0, 
                 4.0 * std::max(1.0, stp0)};

    double stp = 0.0;
    double f0, g0;
    std::tie(f0, g0) = phi(stp);

    stp = stp0;
    int nfev = 0;
    task_value task = task_value::start;

    double f, g;
    std::tie(f, g) = phi(stp);

    do{
        dcsrch(f0, g0, stp, f, g, task, opts);
        if (task == task_value::fg) {
            nfev = nfev + 1;
            std::tie(f, g) = phi(stp);
        }
        if (task == task_value::convergence)
            break;
    }while (true);
    std::cout << std::scientific 
              << std::setw(16) << stp0 
              << std::setw(16) << task 
              << std::setw(4) << nfev 
              << std::setw(16) << stp  
              << std::setw(16) << g << std::endl;    
}

int main(int argc, char** argv){

    const double a0 = 1e-03;
    const double factor = 100;
    const int n = 4;
    // line_search(function, a0, mu, eta)
    std::cout << std::endl << "Table 5.1" << std::endl;
    for (int i = 0; i < n; ++i){
        more_thuente_line_search(std::bind(phi51, std::placeholders::_1, 2.0), std::pow(factor, i) * a0, 1e-03, 1e-01);        
    }

    std::cout << std::endl << "Table 5.2" << std::endl;
    for (int i = 0; i < n; ++i){
        more_thuente_line_search(std::bind(phi52, std::placeholders::_1, 0.004), std::pow(factor, i) * a0, 1e-01, 1e-01);
    }

    std::cout << std::endl << "Table 5.3" << std::endl;
    for (int i = 0; i < n; ++i){
        more_thuente_line_search(std::bind(phi53, std::placeholders::_1, 0.01, 39), std::pow(factor, i) * a0, 1e-01, 1e-01);
    }

    std::cout << std::endl << "Table 5.4" << std::endl;
    for (int i = 0; i < n; ++i){
        more_thuente_line_search(std::bind(phi54, std::placeholders::_1, 0.001, 0.001), std::pow(factor, i) * a0, 1e-03, 1e-03);
    }

    std::cout << std::endl << "Table 5.5" << std::endl;
    for (int i = 0; i < n; ++i){
        more_thuente_line_search(std::bind(phi54, std::placeholders::_1, 0.01, 0.001), std::pow(factor, i) * a0, 1e-03, 1e-03);
    }

    std::cout << std::endl << "Table 5.6" << std::endl;
    for (int i = 0; i < n; ++i){
        more_thuente_line_search(std::bind(phi54, std::placeholders::_1, 0.001, 0.01), std::pow(factor, i) * a0, 1e-03, 1e-03);
    }
    return 0;
}