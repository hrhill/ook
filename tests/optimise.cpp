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

#include "ook/steepest_descent.h"
#include "ook/fletcher_reeves.h"
#include "ook/bfgs.h"
#include "ook/newton.h"
#include "ook/lbfgs.h"
#include "ook/options.h"
#include "ook/stream_observer.h"

#include "ook/test_functions/more_garbow_hillstrom.h"

namespace ublas = boost::numeric::ublas;

typedef ublas::vector<double> vector_type;
typedef ublas::matrix<double, ublas::column_major> matrix_type;

template <typename V, typename M>
using test_function_types = boost::mpl::list<
ook::test_functions::rosenbrock<V, M>,
ook::test_functions::freudenstein_roth<V, M>//,
//ook::test_functions::powell_badly_scaled<V, M>
>;

typedef test_function_types<vector_type, matrix_type> ublas_function_types;

template <typename F, typename X>
struct
gradient_only_wrapper{
    gradient_only_wrapper(F f)
    :
        func(f)
    {}

    std::tuple<double, X>
    operator()(const X& x) const {
        const int n = x.size();
        double f;
        X df(n);
        matrix_type d2f(n, n);

        std::tie(f, df, d2f) = func(x);
        return std::make_tuple(f, df);
    }

    F func;
};

template <typename Function, typename Optimiser>
void
run_gradient_based_optimiser(Function, Optimiser optimiser)
{
    const double epsilon = std::numeric_limits<double>::epsilon();
    ook::options opts{1e-03, 9e-01, epsilon, 0.0, 4.0};

    typedef Function test_function;
    test_function objective_function;

    vector_type x(test_function::n, 0.0);
    std::copy(test_function::x0.begin(), test_function::x0.end(), x.begin());

    vector_type minima(test_function::n, 0.0);
    std::copy(test_function::minima.begin(), test_function::minima.end(), minima.begin());

    gradient_only_wrapper<Function, vector_type> wrapper(objective_function);
    ook::stream_observer<std::ostream> obs(std::cout);
    auto soln = optimiser(wrapper, x, opts, obs);
    BOOST_CHECK_EQUAL(std::get<0>(soln), ook::message::convergence);

    // Evaluate function at minima, check proximity
    /*auto x_min = std::get<1>(soln);
    double f_min;
    vector_type df;
    std::tie(f_min, df) = wrapper(x_min);
    BOOST_CHECK(fabs(f_min - test_function::f_min) <=  1e-08);
    BOOST_CHECK(ook::norm_infinity(x_min - minima) <= 1e-04);
    */
}

template <typename Function>
int
test_gradient_based_optimisers()
{
    std::cout << "steepest_descent" << std::endl;
    run_gradient_based_optimiser(Function(),
            ook::steepest_descent<gradient_only_wrapper<Function, vector_type>,
                                vector_type,
                                ook::options,
                                ook::stream_observer<std::ostream>>);

    std::cout << "fletcher_reeves" << std::endl;
    run_gradient_based_optimiser(Function(),
            ook::fletcher_reeves<gradient_only_wrapper<Function, vector_type>,
                                vector_type,
                                ook::options,
                                ook::stream_observer<std::ostream>>);

    std::cout << "lbfgs" << std::endl;
    run_gradient_based_optimiser(Function(),
            ook::lbfgs<gradient_only_wrapper<Function, vector_type>,
                                vector_type,
                                ook::options,
                                ook::stream_observer<std::ostream>>);

    std::cout << "bfgs" << std::endl;
    run_gradient_based_optimiser(Function(),
            ook::bfgs<gradient_only_wrapper<Function, vector_type>,
                                vector_type,
                                ook::options,
                                ook::stream_observer<std::ostream>>);
    return 0;
}

template <typename Function, typename Optimiser>
void
run_hessian_based_optimiser(Function, Optimiser optimiser)
{
    const double epsilon = std::numeric_limits<double>::epsilon();
    ook::options opts{1e-03, 9e-01, epsilon, 0.0, 4.0};

    Function objective_function;

    vector_type x(Function::n, 0.0);
    std::copy(Function::x0.begin(), Function::x0.end(), x.begin());

    vector_type minima(Function::n, 0.0);
    std::copy(Function::local_minima.begin(), Function::local_minima.end(), minima.begin());

    ook::stream_observer<std::ostream> obs(std::cout);
    auto soln = optimiser(objective_function, x, opts, obs);
    BOOST_CHECK_EQUAL(std::get<0>(soln), ook::message::convergence);
    // Evaluate function at minima, check proximity
    auto x_min = std::get<1>(soln);
    double f_min;
    vector_type df(Function::n);
    matrix_type d2f(Function::n, Function::n);

    std::tie(f_min, df, d2f) = objective_function(x_min);
    //BOOST_CHECK(fabs(f_min - test_function::f_min) <= 1e-08);
    BOOST_CHECK(ook::norm_infinity(x_min - minima) <= 1e-04);
}

template <typename Function>
int test_hessian_based_optimisers()
{
    typedef typename Function::vector_type vector_type;
    std::cout << "newton" << std::endl;
    run_hessian_based_optimiser(Function(),
            ook::newton<Function, vector_type, ook::options, ook::stream_observer<std::ostream>>);
    return 0;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(gradient_based_optimisers, T, ublas_function_types){
    BOOST_CHECK_EQUAL(test_gradient_based_optimisers<T>(), 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(hessian_based_optimisers, T, ublas_function_types){
    BOOST_CHECK_EQUAL(test_hessian_based_optimisers<T>(), 0);
}
