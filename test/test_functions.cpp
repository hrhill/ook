#define BOOST_TEST_MODULE test functions
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <boost/mpl/list.hpp>

#include <iostream>

#include "ook/vector.hpp"
#include "ook/matrix.hpp"
#include "ook/test_functions/parabola.hpp"
#include "ook/test_functions/more_garbow_hillstrom.hpp"

using namespace ook::test_functions;

using test_function_types = boost::mpl::list<
parabola,
rosenbrock,
freudenstein_roth/*,
powell_badly_scaled */
>;

template <typename Function>
int
test_function_specification()
{
    typedef Function test_function;

    test_function objective_function;

    ook::vector x0(test_function::n);
    ook::vector minima(test_function::n);

    std::copy(test_function::x0.begin(), test_function::x0.end(), x0.begin());
    std::copy(test_function::minima.begin(), test_function::minima.end(), minima.begin());

    double f_min;
    ook::vector df(test_function::n);
    ook::matrix d2f(test_function::n, test_function::n);
    std::tie(f_min, df, d2f) = objective_function(minima);

    BOOST_CHECK(abs(f_min - test_function::f_min) <= test_function::tolerance);
    BOOST_CHECK(ook::norm_inf(df) <= test_function::tolerance);
    return 0;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_functions_test, T, test_function_types)
{
    BOOST_CHECK_EQUAL(test_function_specification<T>(), 0);
}
