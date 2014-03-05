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
#include "ook/line_search/backtracking.h"


typedef boost::mpl::list<
    ook::line_search::more_thuente,
    ook::line_search::backtracking
> line_search_types;

std::tuple<double, double>
quadratic(double x)
{
    return std::make_tuple(4 * std::pow(x - 0.5, 2),
                           8 * (x - 0.5));
}

template <typename LineSearch>
int
test_quadratic()
{
    ook::options<double> opts;
    double phi0, dphi0;
    std::tie(phi0, dphi0) = quadratic(0.0);

    double a = 1.0;
    double phia, dphia;
    ook::message msg;
    std::tie(msg, a, phia, dphia) = LineSearch::search(quadratic, phi0, dphi0, a, opts);
    BOOST_CHECK_EQUAL(msg, ook::message::convergence);
    BOOST_CHECK(a < 1.0 && a > 0);
    return 0;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(line_search, T, line_search_types){
    BOOST_CHECK_EQUAL(test_quadratic<T>(), 0);
}
