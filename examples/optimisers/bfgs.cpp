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

#include "norms.h"
#include "options.h"
#include "state_value.h"
#include "line_search/more_thuente.h"
#include "test_functions/more_garbow_hillstrom/rosenbrock.h"

#include "factorisations/gmw81.h"
#include "line_search_methods/steepest_descent.h"
#include "line_search_methods/fletcher_reeves.h"
#include "line_search_methods/bfgs.h"
#include "line_search_methods/newton.h"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> vector_t;
typedef ublas::matrix<double, ublas::column_major> matrix_t;

template <typename Matrix>
Matrix
convert_to_cholesky(const Matrix& LD)
{
    int n = LD.size1();
    Matrix L(LD);
    Matrix D(n, n, 0);

    for (int i = 0; i < n; ++i){
        D(i, i) = sqrt(L(i, i));
        L(i, i) = 1.0;
    }
    return boost::numeric::ublas::prod(L, D);
}

template <typename Matrix, typename Vector>
Vector
solve(Matrix A, const Vector& b)
{
    Matrix LD = ook::factorisations::gmw81(A);
    Matrix L = convert_to_cholesky(LD);

    ublas::symmetric_adaptor<Matrix, ublas::lower> sa(L);    
    Matrix b1(b.size(), 1);

    boost::numeric::ublas::column(b1, 0) = b;
    boost::numeric::bindings::lapack::potrs(sa, b1);

    return ublas::column(b1, 0);
}

template <typename F, typename X, typename Options>
std::tuple<ook::state_value, X>
newton(F objective_function, const X& x0, const Options& opts)
{
    typedef typename X::value_type real_type;

    X x(x0);
    X dfx;
    real_type fx, dfx_dot_p;
    matrix_t d2fx(x0.size(), x0.size());
    // Evaluate at initial point
    std::tie(fx, dfx, d2fx) = objective_function(x0);

    uint iteration = 0;
    uint nfev_total = 0;
    ook::state_value value;

    do {
        // Choose descent direction
        X p = -solve(d2fx, dfx);

        real_type a = 1.0;
        uint nfev = 0;
        dfx_dot_p = std::inner_product(dfx.begin(), dfx.end(), p.begin(), real_type(0.0)); 
        // do line search
        // take a reference to the state variable, ensuring that fx and dfx get updated
        // the line search call will take a fresh copy
        auto phi = [&nfev, &fx, &dfx, &d2fx, &dfx_dot_p, x, p, objective_function](const real_type& a){
            ++nfev;
            std::tie(fx, dfx, d2fx) = objective_function(static_cast<const X&>(x + a * p));
            dfx_dot_p = std::inner_product(dfx.begin(), dfx.end(), p.begin(), real_type(0.0));
            return std::make_pair(fx, dfx_dot_p);
        };

        std::tie(value, a) = ook::line_search::more_thuente(phi, fx, dfx_dot_p, a, opts);

        X dx(a * p);
        x += dx;
        nfev_total += nfev;
        ++iteration;

        ook::detail::report(iteration, nfev_total, nfev, a, fx, dfx, dx);

        if (ook::norm_infinity(dfx) < 1e-08){
            value = ook::state_value::convergence;
            ook::detail::final_report(nfev_total, fx, dfx, dx);                        
            break;
        }

    } while(true);
 
    return std::make_pair(value, x);
}

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
    ook::options opts{1e-03, 9e-01, epsilon, 0.0, 4.0 * std::max(1.0, 1e-03)};

    typedef ook::test_functions::rosenbrock<vector_t, matrix_t> test_function;
    test_function objective_function;
    gradient_only_wrapper<test_function, vector_t> wrapper(objective_function);

    vector_t x(test_function::n, 0.0);
    std::copy(test_function::x0.begin(), test_function::x0.end(), x.begin());

    std::cout << "steepest_descent\n";
    auto soln = ook::steepest_descent(wrapper, x, opts);
    std::cout << std::get<0>(soln) << "\n" << std::get<1>(soln) << std::endl;

    std::cout << "fletcher_reeves\n";
    soln = ook::fletcher_reeves(wrapper, x, opts);
    std::cout << std::get<0>(soln) << "\n" << std::get<1>(soln) << std::endl;

    std::cout << "bfgs\n";
    soln = ook::bfgs(wrapper, x, opts);
    std::cout << std::get<0>(soln) << "\n" << std::get<1>(soln) << std::endl;

    std::cout << "newton\n";
    soln = ook::newton(objective_function, x, opts);
    std::cout << std::get<0>(soln) << "\n" << std::get<1>(soln) << std::endl;
}

