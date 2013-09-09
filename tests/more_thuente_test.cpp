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

#include "more_thuente.h"
#include "line_search_functions.h"

int main(int argc, char** argv){

    const double a0 = 1e-03;
    const double factor = 100;
    const int n = 4;
    const double epsilon = std::numeric_limits<double>::epsilon();

    std::cout << std::endl << "Table 5.1" << std::endl;
    for (int i = 0; i < n; ++i){
        const double stp0 = std::pow(factor, i) * a0;
        options opts{1e-03, 1e-01, epsilon, 0.0, 4.0 * std::max(1.0, stp0)};
        auto s = more_thuente_line_search(std::bind(phi51, std::placeholders::_1, 2.0), stp0, opts);        
        std::cout << std::scientific 
              << std::setw(16) << stp0
              << std::setw(16) << s.task
              << std::setw(4) << s.nfev 
              << std::setw(16) << s.stp  
              << std::setw(16) << s.g << std::endl;    

    }

    std::cout << std::endl << "Table 5.2" << std::endl;
    for (int i = 0; i < n; ++i){
        const double stp0 = std::pow(factor, i) * a0;
        options opts{1e-01, 1e-01, epsilon, 0.0, 4.0 * std::max(1.0, stp0)};        
        auto s =more_thuente_line_search(std::bind(phi52, std::placeholders::_1, 0.004), stp0, opts);
        std::cout << std::scientific 
              << std::setw(16) << stp0
              << std::setw(16) << s.task
              << std::setw(4) << s.nfev 
              << std::setw(16) << s.stp  
              << std::setw(16) << s.g << std::endl;          
    }

    std::cout << std::endl << "Table 5.3" << std::endl;
    for (int i = 0; i < n; ++i){
        const double stp0 = std::pow(factor, i) * a0;
        options opts{1e-01, 1e-01, epsilon, 0.0, 4.0 * std::max(1.0, stp0)};                
        auto s =more_thuente_line_search(std::bind(phi53, std::placeholders::_1, 0.01, 39), stp0, opts);
        std::cout << std::scientific 
              << std::setw(16) << stp0 
              << std::setw(16) << s.task
              << std::setw(4) << s.nfev 
              << std::setw(16) << s.stp  
              << std::setw(16) << s.g << std::endl;          
    }
    std::cout << std::endl << "Table 5.4" << std::endl;
    for (int i = 0; i < n; ++i){
        const double stp0 = std::pow(factor, i) * a0;
        options opts{1e-03, 1e-03, epsilon, 0.0, 4.0 * std::max(1.0, stp0)};        
        auto s =more_thuente_line_search(std::bind(phi54, std::placeholders::_1, 0.001, 0.001), stp0, opts);
        std::cout << std::scientific 
              << std::setw(16) << stp0 
              << std::setw(16) << s.task
              << std::setw(4) << s.nfev 
              << std::setw(16) << s.stp  
              << std::setw(16) << s.g << std::endl;          
    }

    std::cout << std::endl << "Table 5.5" << std::endl;
    for (int i = 0; i < n; ++i){
        const double stp0 = std::pow(factor, i) * a0;
        options opts{1e-03, 1e-03, epsilon, 0.0, 4.0 * std::max(1.0, stp0)};                
        auto s =more_thuente_line_search(std::bind(phi54, std::placeholders::_1, 0.01, 0.001), stp0, opts);
        std::cout << std::scientific 
              << std::setw(16) << stp0
              << std::setw(16) << s.task
              << std::setw(4) << s.nfev 
              << std::setw(16) << s.stp  
              << std::setw(16) << s.g << std::endl;          
    }

    std::cout << std::endl << "Table 5.6" << std::endl;
    for (int i = 0; i < n; ++i){
        const double stp0 = std::pow(factor, i) * a0;
        options opts{1e-03, 1e-03, epsilon, 0.0, 4.0 * std::max(1.0, stp0)};                        
        auto s =more_thuente_line_search(std::bind(phi54, std::placeholders::_1, 0.001, 0.01), stp0, opts);
        std::cout << std::scientific 
              << std::setw(16) << stp0
              << std::setw(16) << s.task
              << std::setw(4) << s.nfev 
              << std::setw(16) << s.stp  
              << std::setw(16) << s.g << std::endl;          
    }
    return 0;
}
