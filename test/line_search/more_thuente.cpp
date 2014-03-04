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

#include "ook/options.h"
#include "ook/message.h"
#include "ook/line_search/more_thuente.h"
#include "ook/test_functions/line_search.h"


BOOST_AUTO_TEST_CASE(safe_guarded_value)
{
    auto linear = [](const double& x){
        std::cout << x << std::endl;
        return ook::test_functions::linear(x, -1.0, 0.0);
    };
/*
    ook::options<double> opts;
    double phi0, dphi0;
    std::tie(phi0, dphi0) = linear(0.0);

    double a = 1.0;
    double phia, dphia;
    ook::message msg;
    std::tie(msg, a, phia, dphia) = ook::line_search::more_thuente()(linear, phi0, dphi0, a, opts);
    BOOST_CHECK_EQUAL(a, 5.0);
    BOOST_CHECK_EQUAL(phia, -5.0);
    BOOST_CHECK_EQUAL(dphia, -1.0);
*/
}
