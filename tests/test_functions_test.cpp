/// \file optimise_test.cpp
#include <iostream>
#include <string>
#include <limits>
#include <random>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>

#include <boost/mpl/list.hpp>

#define BOOST_TEST_MODULE optimise
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include <boost/mpl/list.hpp>

#include "norms.h"
#include "test_functions/more_garbow_hillstrom.h"

using namespace ook::test_functions;

template <typename T>
using test_function_types = boost::mpl::list<
rosenbrock<T>//,
//freudenstein_roth<T>,
//powell_badly_scaled<T>
>;

typedef test_function_types<boost::numeric::ublas::vector<double>> ublas_function_types;

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
    typedef typename vector_type::value_type real_type;
    typedef Function test_function;

    test_function objective_function;

    vector_type x0(test_function::n);
    vector_type minima(test_function::n);

    std::copy(test_function::x0.begin(), test_function::x0.end(), x0.begin());
    std::copy(test_function::minima.begin(), test_function::minima.end(), minima.begin());    

    real_type f_min;
    vector_type df(test_function::n);
    std::tie(f_min, df) = objective_function(minima);

    BOOST_CHECK_CLOSE(f_min, test_function::f_min, test_function::tolerance);
    BOOST_CHECK_SMALL(ook::norm_infinity(df), test_function::tolerance);    
    return 0;
}

// Zero solution
BOOST_AUTO_TEST_CASE_TEMPLATE(test_functions_test, T, ublas_function_types){
    BOOST_CHECK_EQUAL(test_function_specification<T>(), 0);
}