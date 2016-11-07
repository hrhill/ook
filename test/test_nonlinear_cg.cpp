#define BOOST_TEST_MODULE nonlinear_cg

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>

#include <boost/mpl/list.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include "ook/nonlinear_cg.hpp"

// List of beta update methods
typedef boost::mpl::
    list<ook::beta::fr, ook::beta::pr, ook::beta::hs, ook::beta::dy>
        betas;

struct state
{
    state() : dfx(2) {}
    ook::vector dfx;
};

BOOST_AUTO_TEST_CASE_TEMPLATE(nonlinear_cg_descent_direction, T, betas)
{
    // First step is a descent direction.
    state s;
    s.dfx[0] = 1.234;
    s.dfx[1] = 5.678;
    ook::nonlinear_cg_impl<T> scheme(s);
    auto dd = scheme.descent_direction(s);
    BOOST_CHECK_EQUAL(dd[0], -s.dfx[0]);
    BOOST_CHECK_EQUAL(dd[1], -s.dfx[1]);
}
