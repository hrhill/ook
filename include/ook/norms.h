#ifndef OOK_NORMS_H_
#define OOK_NORMS_H_

#include <algorithm>
#include <cmath>

namespace ook{

/// \brief \f$ l_1 \f$ norm.
template <typename T>
typename T::value_type
norm_1(const T& x){
    T r(0.0);
    for (const auto& xi : x)
        r += fabs(xi);
    return r;
}

/// \brief \f$ l_2  \f$ norm.
template <typename T>
typename T::value_type
norm_2(const T& x){
    T r(0.0);
    for (const auto& xi : x)
        r += xi * xi;
    return sqrt(r);
}

/// \brief \f$ l_{\infty} \f$ norm.
template <typename T>
typename T::value_type
norm_infinity(const T& x){
    T r(0.0);
    for (const auto& xi : x){
        const T fxi = fabs(xi);
        if(fxi > r)
            r = fxi;
    }
    return r;
}

}
#endif
