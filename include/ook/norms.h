// Copyright 2013 Harry Hill
//
// This file is part of ook.
//
// ook is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ook is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with ook.  If not, see <http://www.gnu.org/licenses/>.

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
