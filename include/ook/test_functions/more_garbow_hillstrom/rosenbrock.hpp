#ifndef OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_ROSENBROCK_HPP_
#define OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_ROSENBROCK_HPP_

#include <tuple>
#include <limits>
#include <vector>

namespace ook{
namespace test_functions{

template <typename Vector, typename Matrix>
struct rosenbrock
{
    typedef Vector vector_type;
    typedef Matrix matrix_type;

    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type, matrix_type>
    operator()(const vector_type& x) const
    {
        const double x1 = x(0);
        const double x2 = x(1);
        const double f1 = 10 * (x2 - x1 * x1);
        const double f2 = 1- x1;
        const double f = f1 * f1 + f2 * f2;

        vector_type df(2, 1.0);
        df(0) = - 2 * f1 * 20 * x1 - 2 * f2;
        df(1) = 2 * f1 * 10;

        matrix_type d2f(2, 2, 0.0);
        d2f(0, 0) = 2 - 400 * x2 + 1200 * x1 * x1;
        d2f(1, 0) = d2f(0, 1) = -400 * x1;
        d2f(1, 1) = 200;

        return std::make_tuple(f, df, d2f);
    }

    static const int n = 2;
    static const int m = 2;
    static real_type f_min;
    static real_type tolerance;
    static std::vector<real_type> minima;
    static std::vector<real_type> local_minima;
    static std::vector<real_type> x0;
};

template <typename Vector, typename Matrix>
typename Vector::value_type
rosenbrock<Vector, Matrix>::f_min = 0.0;

template <typename Vector, typename Matrix>
typename Vector::value_type
rosenbrock<Vector, Matrix>::tolerance = std::numeric_limits<typename Vector::value_type>::epsilon();

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type>
rosenbrock<Vector, Matrix>::minima = {1.0, 1.0};

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type>
rosenbrock<Vector, Matrix>::local_minima = {1.0, 1.0};


template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type>
rosenbrock<Vector, Matrix>::x0 = {-1.2, 1.0};


} // ns test_functions
} // ns ook

#endif
