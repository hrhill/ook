#ifndef LINALG_OPERATIONS_BLAZE_HPP_
#define LINALG_OPERATIONS_BLAZE_HPP_

#if HAVE_BLAZE

#include "linalg/std_traits.hpp"

#include <blaze/Math.h>

namespace linalg{

template <
    template <typename, bool> class V,
    typename T,
    bool SO>
std::size_t
size(const V<T, SO>& v)
{
    return v.size();
}

template <
    template <typename, bool> class Matrix,
    typename T,
    bool SO>
size_t
num_rows(const Matrix<T, SO>& m)
{
    return m.rows();
}

template <
    template <typename, bool> class Matrix,
    typename T,
    bool SO>
size_t
num_cols(const Matrix<T, SO>& m)
{
    return m.columns();
}

template <
    template <typename, bool> class Matrix,
    typename T,
    bool SO>
auto
column(Matrix<T, SO>& a, size_t idx) -> decltype(blaze::column(a, idx))
{
    return blaze::column(a, idx);
}

template <
    template <typename, bool> class Matrix,
    typename T,
    bool SO>
auto
row(Matrix<T, SO>& a, size_t idx) -> decltype(blaze::row(a, idx))
{
    return blaze::row(a, idx);
}

template <
    template <typename, bool> class Matrix,
    typename T,
    bool SO>
Matrix<T, SO>
trans(const Matrix<T, SO>& M)
{
    return blaze::trans(M);
}

template<
    template <typename, bool> class Matrix1,
    typename T1, bool SO1,
    template <typename, bool> class Matrix2,
    typename T2, bool SO2>
void
swap(blaze::DenseRow<Matrix1<T1, SO1>>& r1, blaze::DenseRow<Matrix2<T2, SO2>>& r2)
{
    const int n = blaze::size(r1);
    std::vector<T1> tmp(n);
    for (int i = 0; i < n; ++i)
        tmp[i] = r1[i];
    for (int i = 0; i < n; ++i)
        r1[i] = r2[i];
    for (int i = 0; i < n; ++i)
        r2[i] = tmp[i];
}

/// \brief Inner product
template <typename T, bool TF1, bool TF2>
auto inner_prod(const blaze::DynamicVector<T, TF1>& x, const blaze::DynamicVector<T, TF2>& y)
{
    typedef remove_const_reference<decltype(x[0])> value_type;
    assert(linalg::size(x) == linalg::size(y));
    return std::inner_product(
                x.begin(), x.end(), y.begin(), value_type(0.0));
}

/// \brief \f$ l_1 \f$ norm.
template <typename T, bool TF>
auto norm_1(const blaze::DynamicVector<T, TF>& x)
{
    typedef remove_const_reference<decltype(x[0])> value_type;
    return std::accumulate(x.begin(), x.end(), value_type(0.0),
            [](const value_type& init, const value_type& xi){
                return init + fabs(xi);
            });
}

/// \brief \f$ l_2  \f$ norm.
template <typename T, bool TF>
auto norm_2(const blaze::DynamicVector<T, TF>& x)
{
    return sqrt(inner_prod(x, x));
}

/// \brief \f$ l_p  \f$ norm.
template <typename T, bool TF>
auto norm_p(const blaze::DynamicVector<T, TF>& x, int p)
{
    assert(p > 0);
    typedef remove_const_reference<decltype(x[0])> value_type;
    value_type r(0.0);
    for (const auto& xi : x)
        r += std::pow(xi, p);
    return exp(log(r)/p);
}

/// \brief \f$ l_{\infty} \f$ norm.
template <typename T, bool TF>
auto norm_infinity(const blaze::DynamicVector<T, TF>& x)
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

#endif
