#ifndef OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_POWELL_BADLY_SCALED_H_
#define OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_POWELL_BADLY_SCALED_H_

#include <tuple>
#include <limits>
#include <vector>

namespace ook{
namespace test_functions{

template <typename Vector, typename Matrix>
struct powell_badly_scaled
{
    typedef Vector vector_type;
    typedef Matrix matrix_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type, matrix_type>
    operator()(const vector_type& x) const
    {
        const real_type x1 = x(0);
        const real_type x2 = x(1);
        const real_type e1 = exp(-x1);
        const real_type e2 = exp(-x2);

        const real_type a = 1e04;
        const real_type f1 = a * x1 * x2 - 1;
        const real_type f2 = e1 + e2 - 1.0001;

        const real_type f = f1 * f1 + f2 * f2;
        vector_type df(2, 0.0);
        matrix_type d2f(2, 2, 0.0);

        const real_type df1x1 = a * x2;
        const real_type df1x2 = a * x1;
        const real_type df2x1 = -e1;
        const real_type df2x2 = -e2;

        const real_type d2f1x1 = 0;
        const real_type d2f1x2 = 0;
        const real_type d2f1x1x2 = a;
        const real_type d2f2x1 = e1;
        const real_type d2f2x2 = e2;
        const real_type d2f2x1x2 = 0;
        df(0) = 2 * f1 * df1x1 + 2 * f2 * df2x1;
        df(1) = 2 * f1 * df1x2 + 2 * f2 * df2x2;

        d2f(0, 0) = 2 * df1x1 * df1x1 + 2 * f1 * d2f1x1 + 2 * df2x1 * df2x1 + 2 * f2 * d2f2x1;
        d2f(1, 1) = 2 * df1x2 * df1x2 + 2 * f1 * d2f1x2 + 2 * df2x2 * df2x2 + 2 * f2 * d2f2x2;
        d2f(0, 1) = d2f(1, 0) = 2 * df1x1 * df1x2 + + 2 * f1 * d2f1x1x2 + 2 * df2x2 * df2x2 + 2 * f2 * d2f2x1x2;

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
powell_badly_scaled<Vector, Matrix>::f_min = 0.0;

template <typename Vector, typename Matrix>
typename Vector::value_type
powell_badly_scaled<Vector, Matrix>::tolerance = std::numeric_limits<typename Vector::value_type>::epsilon();

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type>
powell_badly_scaled<Vector, Matrix>::minima = {1.098159e-05,9.106147};

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type>
powell_badly_scaled<Vector, Matrix>::local_minima = {1.098159e-05,9.106147};

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type>
powell_badly_scaled<Vector, Matrix>::x0 = {0.0,  1.0};


} // ns test_functions
} // ns ook

#endif
