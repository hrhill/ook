/// \file more_thuente_test.cpp
#include <iostream>
#include <string>
#include <limits>
#include <random>
#include <tuple>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>

#include <boost/mpl/list.hpp>

#define BOOST_TEST_MODULE more_thuente
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include "state.h"
#include "options.h"
#include "line_search/more_thuente.h"
#include "line_search_functions.h"

template <typename ObjectiveFunction, typename Options>
void
do_search(ObjectiveFunction obj, const double stp0, const Options& opts)
{
    int nfev = 0;
    double phi0, dphi0, phix, dphix;

    auto phi = [&nfev, &phix, &dphix, obj](const double x){
                        ++nfev;
                        std::tie(phix, dphix) = obj(x);
                        return std::make_tuple(phix, dphix);
                    };

    std::tie(phi0, dphi0) = obj(0.0);    

    auto soln = ook::line_search::more_thuente(phi, phi0, dphi0, stp0, opts);                    
    std::cout << std::scientific 
              << std::setw(16) << stp0 
              << std::setw(16) << std::get<0>(soln)
              << std::setw(4) << nfev 
              << std::setw(16) << std::get<1>(soln)
              << std::setw(16) << dphix << std::endl;       
}

int main(int argc, char** argv){

    const double a0 = 1e-03;
    const double factor = 100;
    const int n = 4;
    const double epsilon = std::numeric_limits<double>::epsilon();

    std::cout << std::endl << "Table 5.1" << std::endl;
    for (int i = 0; i < n; ++i){
        const double stp0 = std::pow(factor, i) * a0;
        ook::options opts{1e-03, 1e-01, epsilon, 0.0, 4.0 * std::max(1.0, stp0)};
        auto objective_function = std::bind(phi51, std::placeholders::_1, 2.0);
        do_search(objective_function, stp0, opts);
    }

    std::cout << std::endl << "Table 5.2" << std::endl;
    for (int i = 0; i < n; ++i){
        const double stp0 = std::pow(factor, i) * a0;
        ook::options opts{1e-01, 1e-01, epsilon, 0.0, 4.0 * std::max(1.0, stp0)};        
        auto objective_function = std::bind(phi52, std::placeholders::_1, 0.004);
        do_search(objective_function, stp0, opts);
    }

    std::cout << std::endl << "Table 5.3" << std::endl;
    for (int i = 0; i < n; ++i){
        const double stp0 = std::pow(factor, i) * a0;
        ook::options opts{1e-01, 1e-01, epsilon, 0.0, 4.0 * std::max(1.0, stp0)};                
        auto objective_function = std::bind(phi53, std::placeholders::_1, 0.01, 39);
        do_search(objective_function, stp0, opts);     
    }
    std::cout << std::endl << "Table 5.4" << std::endl;
    for (int i = 0; i < n; ++i){
        const double stp0 = std::pow(factor, i) * a0;
        ook::options opts{1e-03, 1e-03, epsilon, 0.0, 4.0 * std::max(1.0, stp0)};        
        auto objective_function = std::bind(phi54, std::placeholders::_1, 0.001, 0.001);
        do_search(objective_function, stp0, opts);        
    }

    std::cout << std::endl << "Table 5.5" << std::endl;
    for (int i = 0; i < n; ++i){
        const double stp0 = std::pow(factor, i) * a0;
        ook::options opts{1e-03, 1e-03, epsilon, 0.0, 4.0 * std::max(1.0, stp0)};                
        auto objective_function = std::bind(phi54, std::placeholders::_1, 0.01, 0.001);
        do_search(objective_function, stp0, opts);        
    }

    std::cout << std::endl << "Table 5.6" << std::endl;
    for (int i = 0; i < n; ++i){
        const double stp0 = std::pow(factor, i) * a0;
        ook::options opts{1e-03, 1e-03, epsilon, 0.0, 4.0 * std::max(1.0, stp0)};                        
        auto objective_function = std::bind(phi54, std::placeholders::_1, 0.001, 0.01);
        do_search(objective_function, stp0, opts);
    }
    return 0;
}
