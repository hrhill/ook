#define BOOST_TEST_MODULE steepest_descent

#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include "ook/steepest_descent.hpp"

#include <boost/numeric/ublas/vector.hpp>

typedef boost::numeric::ublas::vector<double> vector_type;

struct state{
    vector_type dfx;
};

BOOST_AUTO_TEST_CASE(steepest_descent_descent_direction)
{
    ook::steepest_descent_impl<vector_type> scheme(1.0);
    state s;
    s.dfx = vector_type(1, 1.234);
    vector_type dd = scheme.descent_direction(s);
    BOOST_CHECK_EQUAL(dd[0], -1.234);
}
