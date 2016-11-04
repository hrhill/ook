#define BOOST_TEST_MODULE bfgs

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>

#include "ook/bfgs.hpp"

struct state
{
    state() : dfx(2), dx(2) {}
    ook::vector dfx;
    ook::vector dx;
};

BOOST_AUTO_TEST_CASE(bfgs_descent_direction)
{
    // First step is a descent direction.
    state s;
    s.dfx[0] = 1.234;
    s.dfx[1] = 5.678;
    ook::bfgs_impl scheme(s);
    auto dd = scheme.descent_direction(s);
    BOOST_CHECK_EQUAL(dd[0], -s.dfx[0]);
    BOOST_CHECK_EQUAL(dd[1], -s.dfx[1]);
}

BOOST_AUTO_TEST_CASE(bfgs_update)
{
    // First step is a descent direction.
    state s;
    s.dfx[0] = 1;
    s.dfx[1] = 1;
    s.dx = s.dfx;
    ook::bfgs_impl scheme(s);

    // Choose dfx and dx so that y = s = (1, 1)
    s.dfx[0] = 2.0;
    s.dfx[1] = 2.0;
    scheme.update(s);

    // New B matrix should be id
    auto dd = scheme.descent_direction(s);
    BOOST_CHECK_EQUAL(dd[0], -s.dfx[0]);
    BOOST_CHECK_EQUAL(dd[1], -s.dfx[1]);
}
