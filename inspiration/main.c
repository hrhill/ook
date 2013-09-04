#include <iostream>
#include <string>
#include <algorithm>
#include <tuple>
#include <iomanip>

#include <boost/math/constants/constants.hpp>

#include "line_search.h"

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
line_search_test(F phi, const double stp0, const double& mu, const double eta){

    const double factor = 1.0;
    const int ntries = 100;

    double stp = 0.0;
    double f, g;

    std::tie(f, g) = phi(stp);

    double ftol = mu;
    double gtol = eta;
    double xtol = std::numeric_limits<double>::epsilon();

    double stpmin = 0.0;
    double stpmax = 4.0 * std::max(1.0, stp0);

    stp = factor * stp0;
    int nfev = 0;
    std::string task("START");

    do{
        dcsrch_(stp, f, g, ftol, gtol, xtol, task, stpmin, stpmax);
        if (task.find("FG") != std::string::npos) {
            nfev = nfev + 1;
            std::tie(f, g) = phi(stp);
        }
    }while (task.find("CONV") == std::string::npos);
    std::cout << std::scientific 
              << std::setw(16) << stp0 
              << std::setw(16) << task 
              << std::setw(4) << nfev 
              << std::setw(16) << stp  
              << std::setw(16) << g << std::endl;    
}

int main(int argc, char** argv){

    // line_search(function, a0, mu, eta)
    std::cout << std::endl << "Table 5.1" << std::endl;
    line_search_test(std::bind(phi51, std::placeholders::_1, 2.0), 1e-03, 1e-03, 1e-01);
    line_search_test(std::bind(phi51, std::placeholders::_1, 2.0), 1e-01, 1e-03, 1e-01);    
    line_search_test(std::bind(phi51, std::placeholders::_1, 2.0), 1e01, 1e-03, 1e-01);    
    line_search_test(std::bind(phi51, std::placeholders::_1, 2.0), 1e03, 1e-03, 1e-01);    

    std::cout << std::endl << "Table 5.2" << std::endl;
    line_search_test(std::bind(phi52, std::placeholders::_1, 0.004), 1e-03, 1e-01, 1e-01);
    line_search_test(std::bind(phi52, std::placeholders::_1, 0.004), 1e-01, 1e-01, 1e-01);    
    line_search_test(std::bind(phi52, std::placeholders::_1, 0.004), 1e01, 1e-01, 1e-01);
    line_search_test(std::bind(phi52, std::placeholders::_1, 0.004), 1e03, 1e-01, 1e-01);    

    std::cout << "Table 5.3" << std::endl;
    line_search_test(std::bind(phi53, std::placeholders::_1, 0.01, 39), 1e-03, 1e-01, 1e-01);
    line_search_test(std::bind(phi53, std::placeholders::_1, 0.01, 39), 1e-01, 1e-01, 1e-01);    
    line_search_test(std::bind(phi53, std::placeholders::_1, 0.01, 39), 1e01, 1e-01, 1e-01);
    line_search_test(std::bind(phi53, std::placeholders::_1, 0.01, 39), 1e03, 1e-01, 1e-01);    

    std::cout << std::endl << "Table 5.4" << std::endl;
    line_search_test(std::bind(phi54, std::placeholders::_1, 0.001, 0.001), 1e-03, 1e-01, 1e-01);
    line_search_test(std::bind(phi54, std::placeholders::_1, 0.001, 0.001), 1e-01, 1e-01, 1e-01);    
    line_search_test(std::bind(phi54, std::placeholders::_1, 0.001, 0.001), 1e01, 1e-01, 1e-01);
    line_search_test(std::bind(phi54, std::placeholders::_1, 0.001, 0.001), 1e03, 1e-01, 1e-01);    

    std::cout << std::endl << "Table 5.5" << std::endl;
    line_search_test(std::bind(phi54, std::placeholders::_1, 0.01, 0.001), 1e-03, 1e-03, 1e-03);
    line_search_test(std::bind(phi54, std::placeholders::_1, 0.01, 0.001), 1e-01, 1e-03, 1e-03);    
    line_search_test(std::bind(phi54, std::placeholders::_1, 0.01, 0.001), 1e01, 1e-03, 1e-03);
    line_search_test(std::bind(phi54, std::placeholders::_1, 0.01, 0.001), 1e03, 1e-03, 1e-03);    

    std::cout << std::endl << "Table 5.6" << std::endl;
    line_search_test(std::bind(phi54, std::placeholders::_1, 0.001, 0.01), 1e-03, 1e-03, 1e-03);
    line_search_test(std::bind(phi54, std::placeholders::_1, 0.001, 0.01), 1e-01, 1e-03, 1e-03);    
    line_search_test(std::bind(phi54, std::placeholders::_1, 0.001, 0.01), 1e01, 1e-03, 1e-03);
    line_search_test(std::bind(phi54, std::placeholders::_1, 0.001, 0.01), 1e03, 1e-03, 1e-03);   
    return 0;
}