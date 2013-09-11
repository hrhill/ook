#ifndef IFDEF_H_
#define IFDEF_H_

#include <vector>

template <typename Vector>
struct NAME
{
    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        BODY
    }

    static const int n = N;
    static const int m = M;
    static real_type f_min;
    static real_type tolerance;
    static std::vector<real_type> minima;
    static std::vector<real_type> x0;
};

template <typename Vector>
typename Vector::value_type
NAME<Vector>::f_min = FMIN;

template <typename Vector>
typename Vector::value_type
NAME<Vector>::tolerance = TOL;

template <typename Vector>
std::vector<typename Vector::value_type>
NAME<Vector>::minima = MINIMA;

template <typename Vector>
std::vector<typename Vector::value_type>
NAME<Vector>::x0 = XO;

#endif