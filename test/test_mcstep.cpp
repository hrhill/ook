#define BOOST_TEST_MODULE mcstep

#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>

#include "ook/line_search/mmt/mcstep.hpp"

BOOST_AUTO_TEST_CASE(is_closer_check)
{
    BOOST_CHECK(ook::line_search::is_closer(1, 0, 1));
    BOOST_CHECK(!ook::line_search::is_closer(1, 2, 2));
}

BOOST_AUTO_TEST_CASE(quadratic_step_check)
{
    // (x - 1/2)^2
    BOOST_CHECK_EQUAL(ook::line_search::quadratic_step(0.0, 0.25, -1.0, 1.0, 0.25), 0.5);
}

BOOST_AUTO_TEST_CASE(cubic_step_check)
{
    // (x - 1/2)^2
    BOOST_CHECK_EQUAL(ook::line_search::cubic_step(0.0, 0.25, -1.0, 1.0, 0.25, 1.0), 0.5);
}
