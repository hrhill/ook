/// \file optimise_test.cpp
#include <iostream>
#include <string>
#include <limits>
#include <random>
#include <iomanip>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>

#include <boost/mpl/list.hpp>

#define BOOST_TEST_MODULE optimise
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <boost/mpl/list.hpp>

#include "ook/test_functions/parabola.hpp"
#include "ook/test_functions/more_garbow_hillstrom.hpp"

using namespace ook::test_functions;

template <typename V, typename M>
using test_function_types = boost::mpl::list<
parabola<V, M>,
rosenbrock<V, M>,
freudenstein_roth<V, M>//,
//powell_badly_scaled<V, M>
>;

typedef test_function_types<
    boost::numeric::ublas::vector<double>,
    boost::numeric::ublas::matrix<double>
    > ublas_function_types;

template <typename Vector>
std::ostream& operator<<(std::ostream& out, const std::tuple<typename Vector::value_type, Vector>& t)
{
    return out << "(" << std::get<0>(t) << ", " << std::get<1>(t) << ")";
}

template <typename Function>
int
test_function_specification()
{
    typedef typename Function::vector_type vector_type;
    typedef typename Function::matrix_type matrix_type;
    typedef typename vector_type::value_type real_type;

    typedef Function test_function;

    test_function objective_function;

    vector_type x0(test_function::n);
    vector_type minima(test_function::n);

    std::copy(test_function::x0.begin(), test_function::x0.end(), x0.begin());
    std::copy(test_function::minima.begin(), test_function::minima.end(), minima.begin());

    real_type f_min;
    vector_type df(test_function::n);
    matrix_type d2f(test_function::n, test_function::n);
    std::tie(f_min, df, d2f) = objective_function(minima);

    BOOST_CHECK(abs(f_min - test_function::f_min) <= test_function::tolerance);
    BOOST_CHECK(boost::numeric::ublas::norm_inf(df) <= test_function::tolerance);
    return 0;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(test_functions_test, T, ublas_function_types){
    BOOST_CHECK_EQUAL(test_function_specification<T>(), 0);
}
