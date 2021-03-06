#define BOOST_TEST_MODULE line_search

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <limits>
#include <random>
#include <string>
#include <tuple>

#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/timer.hpp>

#include <boost/mpl/list.hpp>

#include "ook/line_search/backtracking.hpp"
#include "ook/line_search/mcsrch.hpp"
#include "ook/line_search/options.hpp"
#include "ook/message.hpp"

typedef boost::mpl::list<ook::line_search::mcsrch,
                         ook::line_search::backtracking>
    line_search_types;

std::tuple<double, double>
quadratic(double x)
{
    const double t = x - 0.5;
    return std::make_tuple(4 * std::pow(t, 2), 8 * t);
}

template <typename LineSearch>
int
test_quadratic()
{
    ook::line_search::options<double> opts;
    double phi0, dphi0;
    std::tie(phi0, dphi0) = quadratic(0.0);

    double a = 1.0;
    double phia, dphia;
    ook::message msg;
    LineSearch search;
    std::tie(msg, a, phia, dphia) = search(quadratic, phi0, dphi0, a, opts);
    BOOST_CHECK_EQUAL(msg, ook::message::convergence);
    BOOST_CHECK(a < 1.0 && a > 0);
    return 0;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(line_search, T, line_search_types)
{
    BOOST_CHECK_EQUAL(test_quadratic<T>(), 0);
}
