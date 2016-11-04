#define BOOST_TEST_MODULE newton

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>

#include "ook/newton.hpp"

struct state
{
    state() : dfx(2), dx(2), H(2, 2, 0.0) {}

    ook::vector dfx;
    ook::vector dx;
    ook::matrix H;
};

BOOST_AUTO_TEST_CASE(newton_descent_direction)
{
    // If H is id, then
    state s;
    s.dfx[0] = 1.234;
    s.dfx[1] = 5.678;
    s.H(0, 0) = 1.0;
    s.H(1, 1) = 1.0;

    ook::newton_impl scheme(s);
    auto dd = scheme.descent_direction(s);
    BOOST_CHECK_EQUAL(dd[0], -s.dfx[0]);
    BOOST_CHECK_EQUAL(dd[1], -s.dfx[1]);
}
