#define BOOST_TEST_MODULE steepest_descent

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>

#include "ook/steepest_descent.hpp"
#include "ook/vector.hpp"

struct state
{
    ook::vector dfx;
};

BOOST_AUTO_TEST_CASE(steepest_descent_descent_direction)
{
    ook::steepest_descent_impl scheme(1.0);
    state s;
    s.dfx = ook::vector(1, 1.234);
    ook::vector dd = scheme.descent_direction(s);
    BOOST_CHECK_EQUAL(dd[0], -1.234);
}
