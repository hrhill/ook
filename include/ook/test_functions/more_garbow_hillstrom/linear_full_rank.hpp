#ifndef OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_LINEAR_FULL_RANK_H_
#define OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_LINEAR_FULL_RANK_H_

#include <tuple>
#include <limits>
#include <vector>

namespace ook{
namespace test_functions{

template <typename Vector, typename Matrix>
struct linear_full_rank
{
    typedef Vector vector_type;
    typedef Vector matrix_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type, matrix_type>
    operator()(const vector_type& x) const
    {
        real_type f(1.0);
        vector_type df(0, 1.0);
        matrix_type d2f(0, 0, 0.0);
        return std::make_tuple(f, df, d2f);
    }

    static const int n = 0;
    static const int m = 0;
    static real_type f_min;
    static real_type tolerance;
    static std::vector<real_type> minima;
    static std::vector<real_type> x0;
};

template <typename Vector, typename Matrix>
typename Vector::value_type
linear_full_rank<Vector, Matrix>::f_min = 0.0;

template <typename Vector, typename Matrix>
typename Vector::value_type
linear_full_rank<Vector, Matrix>::tolerance = std::numeric_limits<typename Vector::value_type>::epsilon();

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type>
linear_full_rank<Vector, Matrix>::minima = {};

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type>
linear_full_rank<Vector, Matrix>::x0 = {};


} // ns test_functions
} // ns ook

#endif
