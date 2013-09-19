#include <iostream>
#include <limits>
#include <iomanip>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/assignment.hpp>

#include "state.h"
#include "norms.h"
#include "options.h"
#include "state_value.h"
#include "line_search/more_thuente.h"
#include "test_functions/more_garbow_hillstrom/rosenbrock.h"

#include "factorisations/gmw81.h"

namespace ublas = boost::numeric::ublas;
typedef ublas::vector<double> vector_t;
typedef ublas::matrix<double> matrix_t;

/*
template <typename F, typename X>
struct
line_search_function{

    typedef typename X::value_type real_type;

    line_search_function(F f, X& x, X& dfx, X& p, real_type& fx, real_type& dfx_dot_p, uint& nfev)
    :
        function_(f), x_(x), dfx_(dfx), p_(p), fx_(fx), dfx_dot_p_(dfx_dot_p), nfev_(nfev)
    {}

    std::tuple<real_type, real_type>
    operator()(const real_type& a)
    {
        ++nfev_;
        std::tie(fx_, dfx_) = function_(x_ + a * p_);
        dfx_dot_p_ = std::inner_product(dfx_.begin(), dfx_.end(), p_.begin(), real_type(0.0));
        return std::make_tuple(fx_, dfx_dot_p_);
    }

    F function_;
    const X& x_;
    X& dfx_;
    const X& p_;

    real_type& fx_;
    real_type& dfx_dot_p_;

    uint nfev_;
};
*/

template <typename F, typename X, typename Options>
std::tuple<ook::state_value, X>
steepest_descent(F objective_function, const X& x0, const Options& opts)
{
    typedef typename X::value_type real_type;

    X x(x0);
    X dfx;
    real_type fx, dfx_dot_p;

    // Evaluate at initial point
    std::tie(fx, dfx) = objective_function(x0);

    uint iteration = 0;
    uint nfev_total = 0;
    ook::state_value value;

    do {
        // Choose descent direction
        X p = -dfx;        
        real_type a = 1.0;
        uint nfev = 0;
        dfx_dot_p = std::inner_product(dfx.begin(), dfx.end(), p.begin(), real_type(0.0)); 

        // do line search
        // take a reference to the state variable, ensuring that fx and dfx get updated
        // the line search call will take a fresh copy
        auto phi = [&nfev, &fx, &dfx, &dfx_dot_p, x, p, objective_function](const real_type& a){
            ++nfev;
            std::tie(fx, dfx) = objective_function(static_cast<const X&>(x + a * p));
            dfx_dot_p = std::inner_product(dfx.begin(), dfx.end(), p.begin(), real_type(0.0));
            return std::make_pair(fx, dfx_dot_p);
        };

        std::tie(value, a) = ook::line_search::more_thuente(phi, fx, dfx_dot_p, a, opts);

        X dx(a * p);
        x += dx;
        nfev_total += nfev;
        std::cout << std::setw(8) << ++iteration 
                 << std::scientific 
                 << std::setw(16) << nfev_total
                 << std::setw(8) << nfev 
                 << std::setw(16) << a
                 << std::setw(16) << fx
                 << std::setw(16) << ook::norm_infinity(dfx)
                 << std::setw(16) << ook::norm_infinity(dx) << std::endl;  

        if (ook::norm_infinity(dfx) < 1e-08){
            value = ook::state_value::convergence;
            break;
        }

    } while(true);
 
    return std::make_pair(value, x);
}

template <typename F, typename X, typename Options>
std::tuple<ook::state_value, X>
fletcher_reeves(F objective_function, const X& x0, const Options& opts)
{
    typedef typename X::value_type real_type;

    const int n = x0.size();
    X x(x0);
    X dfx;
    real_type fx, dfx_dot_p;

    // Evaluate at initial point
    std::tie(fx, dfx) = objective_function(x0);

    uint iteration = 0;
    uint nfev_total = 0;
    ook::state_value value;
    double beta = 0;
    X p(n, 0);

    do {
        X dfx0(dfx);

        // Choose descent direction
        p = -dfx + beta * p;

        // do line search
        real_type a = 1.0;
        uint nfev = 0;
        dfx_dot_p = std::inner_product(dfx.begin(), dfx.end(), p.begin(), real_type(0.0)); 
        // take a reference to the state variable, ensuring that fx and dfx get updated
        // the line search call will take a fresh copy
        auto phi = [&nfev, &fx, &dfx, &dfx_dot_p, x, p, objective_function](const real_type& a){
            ++nfev;
            std::tie(fx, dfx) = objective_function(static_cast<const X&>(x + a * p));
            dfx_dot_p = std::inner_product(dfx.begin(), dfx.end(), p.begin(), real_type(0.0));
            return std::make_pair(fx, dfx_dot_p);
        };

        std::tie(value, a) = ook::line_search::more_thuente(phi, fx, dfx_dot_p, a, opts);

        X dx(a * p);
        x += dx;
        nfev_total += nfev;
        std::cout << std::setw(8) << ++iteration 
                 << std::scientific 
                 << std::setw(16) << nfev_total
                 << std::setw(8) << nfev 
                 << std::setw(16) << a
                 << std::setw(16) << fx
                 << std::setw(16) << ook::norm_infinity(dfx)
                 << std::setw(16) << ook::norm_infinity(dx) << std::endl;  

        if (ook::norm_infinity(dfx) < 1e-08){
            value = ook::state_value::convergence;
            break;
        }

        // Update
        beta = ublas::inner_prod(dfx, dfx)/ublas::inner_prod(dfx0, dfx0);
    } while(true);
 
    return std::make_pair(value, x);
}

/// main optimisation loop
template <typename F, typename X, typename Options>
std::tuple<ook::state_value, X>
bfgs(F objective_function, const X& x0, const Options& opts)
{
    typedef typename X::value_type real_type;

    const int n = x0.size();
    X x(x0);
    X dfx;
    real_type fx, dfx_dot_p;

    // BFGS
    matrix_t H = ublas::identity_matrix<double>(n);

    // Evaluate at initial point
    std::tie(fx, dfx) = objective_function(x0);

    uint iteration = 0;
    uint nfev_total = 0;
    ook::state_value value;

    do {
        vector_t y(-dfx);

        // Choose descent direction
        X p = - ublas::prod(H, dfx);

        // do line search
        real_type a = 1.0;
        uint nfev = 0;
        dfx_dot_p = std::inner_product(dfx.begin(), dfx.end(), p.begin(), real_type(0.0)); 
        // take a reference to the state variable, ensuring that fx and dfx get updated
        // the line search call will take a fresh copy
        auto phi = [&nfev, &fx, &dfx, &dfx_dot_p, x, p, objective_function](const real_type& a){
            ++nfev;
            std::tie(fx, dfx) = objective_function(static_cast<const X&>(x + a * p));
            dfx_dot_p = std::inner_product(dfx.begin(), dfx.end(), p.begin(), real_type(0.0));
            return std::make_pair(fx, dfx_dot_p);
        };

        std::tie(value, a) = ook::line_search::more_thuente(phi, fx, dfx_dot_p, a, opts);

        X dx(a * p);
        x += dx;
        nfev_total += nfev;
        std::cout << std::setw(8) << ++iteration 
                 << std::scientific 
                 << std::setw(16) << nfev_total
                 << std::setw(8) << nfev 
                 << std::setw(16) << a
                 << std::setw(16) << fx
                 << std::setw(16) << ook::norm_infinity(dfx)
                 << std::setw(16) << ook::norm_infinity(dx) << std::endl;  

        if (ook::norm_infinity(dfx) < 1e-08){
            value = ook::state_value::convergence;
            break;
        }
        
        // Update H
        // s = dx
        y += dfx;
        const real_type rho = 1.0/(inner_prod(y, dx));
        matrix_t Z(ublas::identity_matrix<double>(n) - rho * ublas::outer_prod(dx, y));
        matrix_t ss = rho * ublas::outer_prod(dx, dx);
        if (iteration == 1){
            const double hii = inner_prod(dx, dx);
            for (int i = 0; i < n; ++i){
                H(i, i) = hii;
            }
        }
        H = ublas::prod(Z, matrix_t(ublas::prod(H, ublas::trans(Z)))) + ss;

    } while(true);
 
    return std::make_pair(value, x);
}

template <typename Matrix, typename Vector>
Vector
solve(Matrix A, const Vector& b){
    return b;
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
        X p = solve(d2fx, -dfx);

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
        std::cout << std::setw(8) << ++iteration 
                 << std::scientific 
                 << std::setw(16) << nfev_total
                 << std::setw(8) << nfev 
                 << std::setw(16) << a
                 << std::setw(16) << fx
                 << std::setw(16) << ook::norm_infinity(dfx)
                 << std::setw(16) << ook::norm_infinity(dx) << std::endl;  

        if (ook::norm_infinity(dfx) < 1e-08){
            value = ook::state_value::convergence;
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
/*
    const double epsilon = std::numeric_limits<double>::epsilon();    
    ook::options opts{1e-03, 9e-01, epsilon, 0.0, 4.0 * std::max(1.0, 1e-03)};

    typedef ook::test_functions::rosenbrock<vector_t, matrix_t> test_function;
    test_function objective_function;
    gradient_only_wrapper<test_function, vector_t> wrapper(objective_function);

    vector_t x(test_function::n, 0.0);
    std::copy(test_function::x0.begin(), test_function::x0.end(), x.begin());

    auto soln = steepest_descent(wrapper, x, opts);

    std::cout << std::get<0>(soln) << "\n" << std::get<1>(soln) << std::endl;

    soln = fletcher_reeves(wrapper, x, opts);

    std::cout << std::get<0>(soln) << "\n" << std::get<1>(soln) << std::endl;

    soln = bfgs(wrapper, x, opts);

    std::cout << std::get<0>(soln) << "\n" << std::get<1>(soln) << std::endl;

    soln = newton(objective_function, x, opts);

    std::cout << std::get<0>(soln) << "\n" << std::get<1>(soln) << std::endl;
*/

    matrix_t G(3, 3);
    G <<= 1.0, 1.0,         2.0, 
          1.0, nextafter(1.0, 2), 3.0,
          2.0, 3.0,         1.0;

    std::cout << G << std::endl;
    std::cout << ook::factorisations::gmw81(G) << std::endl;
}

