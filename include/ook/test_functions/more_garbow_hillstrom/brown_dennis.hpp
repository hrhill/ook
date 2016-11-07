#ifndef OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_BROWN_DENNIS_HPP_
#define OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_BROWN_DENNIS_HPP_

#include <limits>
#include <tuple>
#include <vector>

namespace ook
{
namespace test_functions
{

template <typename Vector, typename Matrix>
struct brown_dennis
{
    typedef Vector vector_type;
    typedef Vector matrix_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type, matrix_type>
    operator()(const vector_type& x) const
    {
        real_type f(1.0);
        vector_type df(4, 1.0);
        matrix_type d2f(4, 4, 0.0);
        return std::make_tuple(f, df, d2f);
    }

    static const int n = 4;
    static const int m = 20;
    static real_type f_min;
    static real_type tolerance;
    static std::vector<real_type> minima;
    static std::vector<real_type> x0;
};

template <typename Vector, typename Matrix>
typename Vector::value_type brown_dennis<Vector, Matrix>::f_min = 85822.2;

template <typename Vector, typename Matrix>
typename Vector::value_type brown_dennis<Vector, Matrix>::tolerance =
    std::numeric_limits<typename Vector::value_type>::epsilon();

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type> brown_dennis<Vector, Matrix>::minima =
    {};

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type> brown_dennis<Vector, Matrix>::x0 = {
    25, 5, -5, -1};

} // ns test_functions
} // ns ook

#endif
