#ifndef LINALG_NORMS_HPP_
#define LINALG_NORMS_HPP_

#include <cmath>
#include <numeric>

#include <boost/numeric/ublas/vector_expression.hpp>

#include "operations.hpp"
#include "std_traits.hpp"

namespace linalg{

/// \brief Inner product
template <typename T>
auto inner_prod(const T& x, const T& y)
{
    typedef remove_const_reference<decltype(x[0])> value_type;
    assert(size(x) == size(y));
    return std::inner_product(x.begin(), x.end(), y.begin(), value_type(0.0));
}

/// \brief \f$ l_1 \f$ norm.
template <typename T>
auto norm_1(const boost::numeric::ublas::vector_expression<T>& x)
{
    typedef remove_const_reference<decltype(x()[0])> value_type;
    return std::accumulate(x().begin(), x().end(), value_type(0.0),
            [](const value_type& init, const value_type& xi){
                return init + fabs(xi);
            });
}



/// \brief \f$ l_2  \f$ norm.
template <typename T>
auto norm_2(const boost::numeric::ublas::vector_expression<T>& x)
{
    return sqrt(inner_prod(x(), x()));
}

/// \brief \f$ l_p  \f$ norm.
template <typename T>
auto norm_p(const boost::numeric::ublas::vector_expression<T>& x, int p)
{
    assert(p > 0);
    typedef remove_const_reference<decltype(x[0])> value_type;
    value_type r(0.0);
    for (const auto& xi : x)
        r += std::pow(xi, p);
    return exp(log(r)/p);
}

/// \brief \f$ l_{\infty} \f$ norm.
template <typename T>
auto norm_infinity(const boost::numeric::ublas::vector_expression<T>& x)
{
    typedef remove_const_reference<decltype(x[0])> value_type;
    value_type r(0.0);
    for (const auto& xi : x){
        const value_type fxi = fabs(xi);
        if (fxi > r)
            r = fxi;
    }
    return r;
}

} // ns linalg

#endif
