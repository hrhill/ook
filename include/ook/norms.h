#ifndef OOK_NORMS_H_
#define OOK_NORMS_H_

#include <cmath>
#include <numeric>

namespace ook{

/// \brief Inner product
template <typename T>
typename T::value_type
inner_product(const T& x, const T& y){
    typedef typename T::value_type value_type;
    return std::inner_product(x.begin(), x.end(), y.begin(), value_type(0.0));
}

/// \brief \f$ l_1 \f$ norm.
template <typename T>
typename T::value_type
norm_1(const T& x){
    typename T::value_type r(0.0);
    for (const auto& xi : x)
        r += fabs(xi);
    return r;
}

/// \brief \f$ l_2  \f$ norm.
template <typename T>
typename T::value_type
norm_2(const T& x){
    return sqrt(inner_product(x, x));
}

/// \brief \f$ l_p  \f$ norm.
template <typename T>
typename T::value_type
norm_p(const T& x, int p){
    typename T::value_type r(0.0);
    for (const auto& xi : x)
        r += std::pow(xi, p);
    return exp(log(r)/p);
}

/// \brief \f$ l_{\infty} \f$ norm.
template <typename T>
typename T::value_type
norm_infinity(const T& x){
    typename T::value_type r(0.0);
    for (const auto& xi : x){
        const typename T::value_type fxi = fabs(xi);
        if (fxi > r)
            r = fxi;
    }
    return r;
}

} // ns ook

#endif
