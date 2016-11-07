#define BOOST_TEST_MODULE optimise

#include <iostream>
#include <limits>

#include <boost/mpl/list.hpp>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>

#include "ook/bfgs.hpp"
#include "ook/newton.hpp"
#include "ook/nonlinear_cg.hpp"
#include "ook/options.hpp"
#include "ook/steepest_descent.hpp"
#include "ook/stream_observer.hpp"

#include "ook/test_functions/more_garbow_hillstrom.hpp"
#include "ook/test_functions/parabola.hpp"

using test_function_types =
    boost::mpl::list<ook::test_functions::parabola,
                     ook::test_functions::rosenbrock
                     //, ook::test_functions::freudenstein_roth<V, M>
                     //, ook::test_functions::powell_badly_scaled<V, M>
                     >;

template <typename F>
struct gradient_only_wrapper
{
    explicit gradient_only_wrapper(F f) : func(f) {}

    std::tuple<double, ook::vector>
    operator()(const ook::vector& x) const
    {
        const int n = x.size();
        double f;
        ook::vector df(n);
        ook::matrix d2f(n, n);

        std::tie(f, df, d2f) = func(x);
        return std::make_tuple(f, df);
    }

    F func;
};

template <typename Function, typename Optimiser>
void
run_gradient_based_optimiser(Function, Optimiser optimiser)
{
    ook::options<double> opts;

    Function objective_function;

    ook::vector x(Function::n, 0.0);
    std::copy(Function::x0.begin(), Function::x0.end(), x.begin());

    ook::vector minima(Function::n, 0.0);
    std::copy(Function::minima.begin(), Function::minima.end(), minima.begin());

    gradient_only_wrapper<Function> wrapper(objective_function);
    ook::stream_observer<std::ostream> obs(std::cout);
    auto soln = optimiser(wrapper, x, opts, obs);
    BOOST_CHECK_EQUAL(soln.msg, ook::message::convergence);

    // Evaluate function at minima, check proximity
    auto x_min = soln.x;
    double f_min;
    ook::vector df;
    std::tie(f_min, df) = wrapper(x_min);
    BOOST_CHECK(fabs(f_min - Function::f_min) <= 1e-08);
    BOOST_CHECK(ook::norm_inf(x_min - minima) <= 1e-04);
}

template <typename Function>
int
test_gradient_based_optimisers()
{
    /*
        std::cout << "steepest_descent" << std::endl;
        run_gradient_based_optimiser(Function(),
                ook::steepest_descent<gradient_only_wrapper<Function>,
                                    ook::vector,
                                    ook::options<double>,
                                    ook::stream_observer<std::ostream>>);
    */
    std::cout << "fletcher_reeves" << std::endl;
    run_gradient_based_optimiser(
        Function(),
        ook::fletcher_reeves<gradient_only_wrapper<Function>,
                             ook::options<double>,
                             ook::stream_observer<std::ostream>>);

    std::cout << "bfgs" << std::endl;
    run_gradient_based_optimiser(Function(),
                                 ook::bfgs<gradient_only_wrapper<Function>,
                                           ook::options<double>,
                                           ook::stream_observer<std::ostream>>);
    return 0;
}

template <typename Function, typename Optimiser>
void
run_hessian_based_optimiser(Function, Optimiser optimiser)
{
    ook::vector x(Function::n, 0.0);
    std::copy(Function::x0.begin(), Function::x0.end(), x.begin());

    ook::vector minima(Function::n, 0.0);
    std::copy(Function::local_minima.begin(),
              Function::local_minima.end(),
              minima.begin());

    ook::stream_observer<std::ostream> obs(std::cout);

    ook::options<double> opts;
    Function objective_function;
    auto soln = optimiser(objective_function, x, opts, obs);
    BOOST_CHECK_EQUAL(soln.msg, ook::message::convergence);
    BOOST_CHECK(fabs(soln.fx - Function::f_min) <= 1e-08);
    BOOST_CHECK(ook::norm_inf(soln.x - minima) <= 1e-04);
}

template <typename Function>
int
test_hessian_based_optimisers()
{
    std::cout << "newton" << std::endl;
    run_hessian_based_optimiser(
        Function(),
        ook::newton<Function,
                    ook::options<double>,
                    ook::stream_observer<std::ostream>>);
    return 0;
}

BOOST_AUTO_TEST_CASE_TEMPLATE(gradient_based_optimisers, T, test_function_types)
{
    BOOST_CHECK_EQUAL(test_gradient_based_optimisers<T>(), 0);
}

// BOOST_AUTO_TEST_CASE_TEMPLATE(hessian_based_optimisers, T,
// test_function_types){
//    BOOST_CHECK_EQUAL(test_hessian_based_optimisers<T>(), 0);
//}
