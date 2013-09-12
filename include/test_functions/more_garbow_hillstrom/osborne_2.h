#ifndef OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_OSBORNE_2_H_
#define OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_OSBORNE_2_H_

#include <tuple>
#include <limits>
#include <vector>

namespace ook{
namespace test_functions{

template <typename Vector>
struct osborne_2
{
    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        real_type f(1.0);
        vector_type df(11, 1.0);
        return std::make_pair(f, df);
    }

    static const int n = 11;
    static const int m = 65;
    static real_type f_min;
    static real_type tolerance;
    static std::vector<real_type> minima;
    static std::vector<real_type> x0;
};

template <typename Vector>
typename Vector::value_type
osborne_2<Vector>::f_min = 4.01377e-02;

template <typename Vector>
typename Vector::value_type
osborne_2<Vector>::tolerance = std::numeric_limits<typename Vector::value_type>::epsilon();

template <typename Vector>
std::vector<typename Vector::value_type>
osborne_2<Vector>::minima = {};

template <typename Vector>
std::vector<typename Vector::value_type>
osborne_2<Vector>::x0 = {1.3, 0.65, 0.65, 0.7, 0.6, 3, 5, 7, 2, 4.5, 5.5};


} // ns test_functions
} // ns ook

#endif
