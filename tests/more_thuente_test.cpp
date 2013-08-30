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

std::tuple<double, double>
phi(const double& a)
{
    const double s = 0.1;
    if (a >= 0 && a <= 1){
        return std::make_tuple(0.5 * (1 - s) * a * a - a, 
                                     (1 - s) * a - 1);
    }else{
        return std::make_tuple(0.5 * (s - 1) - s * a,
                                -s);
    }
}

BOOST_AUTO_TEST_CASE(compilation_test){
    double phi0, dphi0;
    std::tie(phi0, dphi0) = phi(0.0);
    BOOST_CHECK(false);
}
