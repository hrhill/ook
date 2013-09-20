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

