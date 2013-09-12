#ifndef OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_GULF_RND_H_
#define OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_GULF_RND_H_

#include <tuple>
#include <limits>
#include <vector>

namespace ook{
namespace test_functions{

template <typename Vector>
struct gulf_rnd
{
    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        real_type f(1.0);
        vector_type df(3, 1.0);
        return std::make_pair(f, df);
    }

    static const int n = 3;
    static const int m = 50;
    static real_type f_min;
    static real_type tolerance;
    static std::vector<real_type> minima;
    static std::vector<real_type> x0;
};

template <typename Vector>
typename Vector::value_type
gulf_rnd<Vector>::f_min = 0.0;

template <typename Vector>
typename Vector::value_type
gulf_rnd<Vector>::tolerance = std::numeric_limits<typename Vector::value_type>::epsilon();

template <typename Vector>
std::vector<typename Vector::value_type>
gulf_rnd<Vector>::minima = {50, 25, 1.5};

template <typename Vector>
std::vector<typename Vector::value_type>
gulf_rnd<Vector>::x0 = {5, 2.5, 0.15};


} // ns test_functions
} // ns ook

#endif
