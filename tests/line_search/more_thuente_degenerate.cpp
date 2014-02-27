/// \file optimise_test.cpp
#include <iostream>
#include <string>
#include <limits>
#include <random>

#include <boost/timer.hpp>

#include <boost/mpl/list.hpp>

#define BOOST_TEST_MODULE more_thuente_degenerate
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include "ook/line_search/more_thuente/more_thuente.h"
#include "ook/options.h"
#include "ook/test_functions/line_search.h"

BOOST_AUTO_TEST_CASE(constant_check){

    const double epsilon = std::numeric_limits<double>::epsilon();
    const double stp0 = 1.0;
    ook::options opts{1e-03, 1e-01, epsilon, 0.0, 4.0 * std::max(1.0, stp0)};

    int nfev = 0;
    double phi0, dphi0, phix, dphix;

    const double a =-1.0;
    const double b = 10.0;

    auto phi = [&nfev, &phix, &dphix, a, b](const double x){
                        ++nfev;
                        std::tie(phix, dphix) = linear(x, a, b);
                        return std::make_tuple(phix, dphix);
                    };

    std::tie(phi0, dphi0) = phi(0.0);
    std::tie(phix, dphix) = phi(1.0);

    auto soln = ook::line_search::more_thuente(phi, phi0, dphi0, stp0, opts);

    std::cout << phix << std::endl;
    std::cout << dphix << std::endl;
    std::cout << std::get<0>(soln) << std::endl;
    std::cout << std::get<1>(soln) << std::endl;
}

