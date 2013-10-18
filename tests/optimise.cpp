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

#include "ook/line_search_methods/steepest_descent.h"
#include "ook/line_search_methods/fletcher_reeves.h"
#include "ook/line_search_methods/bfgs.h"
#include "ook/line_search_methods/newton.h"
#include "ook/line_search_methods/lbfgs.h"
#include "ook/options.h"

#include "ook/test_functions/more_garbow_hillstrom.h"

namespace ublas = boost::numeric::ublas;

typedef ublas::vector<double> vector_t;
typedef ublas::matrix<double, ublas::column_major> matrix_t;

template <typename V, typename M>
using test_function_types = boost::mpl::list<
ook::test_functions::rosenbrock<V, M>,
ook::test_functions::freudenstein_roth<V, M>//,
//ook::test_functions::powell_badly_scaled<V, M>
>;

typedef test_function_types<ublas::vector<double>, ublas::matrix<double>> ublas_function_types;

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
        matrix_t d2f(n, n);

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

    vector_t x(test_function::n, 0.0);
    std::copy(test_function::x0.begin(), test_function::x0.end(), x.begin());

    vector_t minima(test_function::n, 0.0);
    std::copy(test_function::minima.begin(), test_function::minima.end(), minima.begin());

    gradient_only_wrapper<Function, vector_t> wrapper(objective_function);
    auto soln = optimiser(wrapper, x, opts, std::cout);
    BOOST_CHECK_EQUAL(std::get<0>(soln), ook::message::convergence);

    // Evaluate function at minima, check proximity 
    /*auto x_min = std::get<1>(soln);
    double f_min;
    vector_t df;
    std::tie(f_min, df) = wrapper(x_min);
    BOOST_CHECK(fabs(f_min - test_function::f_min) <=  1e-08);
    BOOST_CHECK(ook::norm_infinity(x_min - minima) <= 1e-04);        
    */
}

template <typename Function>
int
test_gradient_based_optimisers()
{
    typedef typename Function::vector_type vector_type;


    std::cout << "steepest_descent" << std::endl;
    run_gradient_based_optimiser(Function(), 
            ook::steepest_descent<gradient_only_wrapper<Function, vector_t>, vector_type, ook::options, std::ostream>);
    std::cout << "fletcher_reeves" << std::endl;    
    run_gradient_based_optimiser(Function(), 
            ook::fletcher_reeves<gradient_only_wrapper<Function, vector_t>, vector_type, ook::options, std::ostream>);
    std::cout << "lbfgs" << std::endl;
    run_gradient_based_optimiser(Function(), 
            ook::lbfgs<gradient_only_wrapper<Function, vector_t>, vector_type, ook::options, std::ostream>);
    std::cout << "bfgs" << std::endl;        
    run_gradient_based_optimiser(Function(), 
            ook::bfgs<gradient_only_wrapper<Function, vector_t>, vector_type, ook::options, std::ostream>);
    return 0;
}

template <typename Function, typename Optimiser>
void
run_hessian_based_optimiser(Function, Optimiser optimiser)
{
    typedef Function test_function;
    typedef typename Function::vector_type vector_type;
    typedef typename Function::matrix_type matrix_type;

    const double epsilon = std::numeric_limits<double>::epsilon();    
    ook::options opts{1e-03, 9e-01, epsilon, 0.0, 4.0};

    test_function objective_function;

    vector_type x(test_function::n, 0.0);
    std::copy(test_function::x0.begin(), test_function::x0.end(), x.begin());

    vector_type minima(test_function::n, 0.0);
    std::copy(test_function::local_minima.begin(), test_function::local_minima.end(), minima.begin());

    auto soln = optimiser(objective_function, x, opts, std::cout);
    BOOST_CHECK_EQUAL(std::get<0>(soln), ook::message::convergence);
    // Evaluate function at minima, check proximity 
    auto x_min = std::get<1>(soln);
    double f_min;
    vector_type df(test_function::n);
    matrix_type d2f(test_function::n, test_function::n);

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
            ook::newton<Function, vector_type, ook::options, std::ostream>);
    return 0;   
}

BOOST_AUTO_TEST_CASE_TEMPLATE(gradient_based_optimisers, T, ublas_function_types){
    BOOST_CHECK_EQUAL(test_gradient_based_optimisers<T>(), 0);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(hessian_based_optimisers, T, ublas_function_types){
    BOOST_CHECK_EQUAL(test_hessian_based_optimisers<T>(), 0);
}