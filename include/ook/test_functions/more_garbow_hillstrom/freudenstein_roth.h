#ifndef OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_FREUDENSTEIN_ROTH_H_
#define OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_FREUDENSTEIN_ROTH_H_

#include <tuple>
#include <limits>
#include <vector>

namespace ook{
namespace test_functions{

template <typename Vector, typename Matrix>
struct freudenstein_roth
{
    typedef Vector vector_type;
    typedef Matrix matrix_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type, matrix_type>
    operator()(const vector_type& x) const
    {
        vector_type df(2, 1.0);
        matrix_type d2f(2, 2, 0.0);

        const real_type x1 = x(0);
        const real_type x2 = x(1);

        const real_type f1 = -13 + x1 + ((5 - x2) * x2 - 2) * x2;
        const real_type f2 = -29 + x1 + ((x2 + 1) * x2 - 14) * x2;
        const real_type f = f1 * f1 + f2 * f2;

        const real_type df1x1 = 1;
        const real_type df1x2 = -3 * x2 * x2 + 10 * x2 - 2;
        const real_type d2f1x1 = 0;
        const real_type d2f1x2 = -6 * x2 + 10;
        const real_type d2f1x1x2 = 0;

        const real_type df2x1 = 1;
        const real_type df2x2 = 3 * x2 * x2 + 2 * x2 - 14;
        const real_type d2f2x1 = 0;        
        const real_type d2f2x2 = 6 * x2 + 2;
        const real_type d2f2x1x2 = 0;

        df(0) = 2 * f1 * df1x1 + 2 * f2 * df2x1;
        df(1) = 2 * f1 * df1x2 + 2 * f2 * df2x2;

        d2f(0, 0) = 2 * f1 * d2f1x1 + 2 * df1x1 * df1x1 + 2 * f2 * d2f2x1 + 2 * df2x1 * df2x1;
        d2f(1, 1) = 2 * f1 * d2f1x2 + 2 * df1x2 * df1x2 + 2 * f2 * d2f2x2 + 2 * df2x2 * df2x2;

        d2f(0, 1) = d2f(1, 0) = 2 * df1x1 * df1x2 + 2 * (f1 + f2) * d2f1x1x2 + 2 * df2x1 * df2x2;

        return std::make_tuple(f, df, d2f);
    }

    static const int n = 2;
    static const int m = 2;
    static real_type f_min;
    static real_type tolerance;
    static std::vector<real_type> minima;
    static std::vector<real_type> x0;
    static std::vector<real_type> local_minima;
};

template <typename Vector, typename Matrix>
typename Vector::value_type
freudenstein_roth<Vector, Matrix>::f_min = 0.0;

template <typename Vector, typename Matrix>
typename Vector::value_type
freudenstein_roth<Vector, Matrix>::tolerance = std::numeric_limits<typename Vector::value_type>::epsilon();

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type>
freudenstein_roth<Vector, Matrix>::minima = {5.0, 4.0};

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type>
freudenstein_roth<Vector, Matrix>::local_minima = {11.41278, -0.8968053};

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type>
freudenstein_roth<Vector, Matrix>::x0 = {0.5,  2.0};


} // ns test_functions
} // ns ook

#endif
