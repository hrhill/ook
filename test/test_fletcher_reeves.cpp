#define BOOST_TEST_MODULE fletcher_reeves

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include "ook/fletcher_reeves.hpp"

#include <boost/numeric/ublas/vector.hpp>

typedef boost::numeric::ublas::vector<double> vector_type;

struct state
{
    state()
    :   dfx(2)
    {}
    vector_type dfx;
};

BOOST_AUTO_TEST_CASE(fletcher_reeves_descent_direction)
{
    // First step is a descent direction.
    state s;
    s.dfx[0] = 1.234;
    s.dfx[1] = 5.678;
    ook::detail::fletcher_reeves<vector_type> scheme(s);
    auto dd = scheme.descent_direction(s);
    BOOST_CHECK_EQUAL(dd[0], -s.dfx[0]);
    BOOST_CHECK_EQUAL(dd[1], -s.dfx[1]);
}

BOOST_AUTO_TEST_CASE(fletcher_reeves_update)
{
    // First step is a descent direction.
    state s;
    s.dfx[0] = 1;
    s.dfx[1] = 1;
    ook::detail::fletcher_reeves<vector_type> scheme(s);

    // First step is descent direction
    // stores p
    auto dd = scheme.descent_direction(s);
    scheme.update(s);

    dd = scheme.descent_direction(s);
    BOOST_CHECK_EQUAL(dd[0], -2.0);
    BOOST_CHECK_EQUAL(dd[1], -2.0);
    scheme.update(s);

    dd = scheme.descent_direction(s);
    BOOST_CHECK_EQUAL(dd[0], -3.0);
    BOOST_CHECK_EQUAL(dd[1], -3.0);
    scheme.update(s);
}
