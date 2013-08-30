#ifndef OOK_NORMS_H_
#define OOK_NORMS_H_

#include <algorithm>
#include <cmath>

namespace ook{

namespace detail{
template <typename T, typename Accumulator, typename Reduce>
typename T::value_type
norm(const T& x, Accumulator acc, Reduce reduce){
    typedef typename T::value_type real_type;
    return reduce(std::accumulate(x.begin(), x.end(), real_type(0.0), acc));
}
}

/// \brief \f$ l_1 \f$ norm.
template <typename T>
typename T::value_type
norm_1(const T& x){
    typedef typename T::value_type real_type;   
    auto acc = [](const real_type s, const real_type xi){
                        return s + sqrt(xi * xi);   
                });
    auto reduce = [](const real_type& x){   
                        return x;   
                    };
    return detail::norm(x, acc, reduce);
}

/// \brief \f$ l_2  \f$ norm.
template <typename T>
typename T::value_type
norm_2(const T& x){
    typedef typename T::value_type real_type;   
    auto acc = [](const real_type s, const real_type xi){
                        return s + xi * xi; 
                });
    return detail::norm(x, acc, sqrt);
}

/// \brief \f$ l_{\infty} \f$ norm.
template <typename T>
typename T::value_type
norm_infinity(const T& x){
    typedef typename T::value_type real_type;
    return *std::max_element(x.begin(), x.end(),
                    [](const real_type& xi, const real_type&
 yi){
                        return fabs(xi) < fabs(yi);
                    });
}

}
#endif
