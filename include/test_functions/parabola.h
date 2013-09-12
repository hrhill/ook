#ifndef OOK_TEST_FUNCTIONS_PARABOLA_H_
#define OOK_TEST_FUNCTIONS_PARABOLA_H_

#include <tuple>
#include <limits>
#include <vector>

#include "../norms.h"

namespace ook{
namespace test_functions{

template <typename Vector>
struct parabola
{
    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        real_type f = 0.5 * std::pow(norm_2(x), 2);
        return std::make_pair(f, x);
    }

    static const int n = 4;
    static const int m = 4;
    static real_type f_min;
    static real_type tolerance;
    static std::vector<real_type> minima;
    static std::vector<real_type> x0;
};

template <typename Vector>
typename Vector::value_type
parabola<Vector>::f_min = 0.0;

template <typename Vector>
typename Vector::value_type
parabola<Vector>::tolerance = std::numeric_limits<typename Vector::value_type>::epsilon();

template <typename Vector>
std::vector<typename Vector::value_type>
parabola<Vector>::minima = {0, 0, 0, 0};

template <typename Vector>
std::vector<typename Vector::value_type>
parabola<Vector>::x0 = {10.0, -4.1513, 1.0, 1000.0};


} // ns test_functions
} // ns ook

#endif
