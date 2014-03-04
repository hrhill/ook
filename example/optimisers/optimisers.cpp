#include <iostream>
#include <limits>
#include <iomanip>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/assignment.hpp>

#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/lapack/computational/potrs.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>

#include "ook/norms.h"
#include "ook/options.h"
#include "ook/message.h"
#include "ook/stream_observer.h"

#include "ook/test_functions/more_garbow_hillstrom/rosenbrock.h"

#include "ook/steepest_descent.h"
#include "ook/fletcher_reeves.h"
#include "ook/bfgs.h"
#include "ook/newton.h"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> vector_t;
typedef ublas::matrix<double, ublas::column_major> matrix_t;

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

int main(){

    const double epsilon = std::numeric_limits<double>::epsilon();
    ook::options<double> opts{1e-03, 9e-01, epsilon, 0.0, 4.0 * std::max(1.0, 1e-03)};

    typedef ook::test_functions::rosenbrock<vector_t, matrix_t> test_function;
    test_function objective_function;
    gradient_only_wrapper<test_function, vector_t> wrapper(objective_function);

    vector_t x(test_function::n, 0.0);
    std::copy(test_function::x0.begin(), test_function::x0.end(), x.begin());
    {
        std::cout << "steepest_descent\n";
        ook::stream_observer<std::ostream> obs(std::cout);
        auto soln = ook::steepest_descent(wrapper, x, opts, obs);
        std::cout << std::get<0>(soln) << "\n" << std::get<1>(soln) << std::endl;
    }

    {
        std::cout << "fletcher_reeves\n";
        ook::stream_observer<std::ostream> obs(std::cout);
        auto soln = ook::fletcher_reeves(wrapper, x, opts, obs);
        std::cout << std::get<0>(soln) << "\n" << std::get<1>(soln) << std::endl;
    }

    {
        std::cout << "bfgs\n";
        ook::stream_observer<std::ostream> obs(std::cout);
        auto soln = ook::bfgs(wrapper, x, opts, obs);
        std::cout << std::get<0>(soln) << "\n" << std::get<1>(soln) << std::endl;
    }

    {
        std::cout << "newton\n";
        ook::stream_observer<std::ostream> obs(std::cout);
        auto soln = ook::newton(objective_function, x, opts, obs);
        std::cout << std::get<0>(soln) << "\n" << std::get<1>(soln) << std::endl;
    }
}

