#define BOOST_TEST_MODULE steepest_descent

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include "ook/steepest_descent.hpp"

struct state{
    double dfx;
};

BOOST_AUTO_TEST_CASE(steepest_descent_descent_direction)
{
    ook::detail::steepest_descent scheme(1);
    state s;
    s.dfx = 1.234;
    BOOST_CHECK_EQUAL(scheme.descent_direction(s), -1.234);
}
