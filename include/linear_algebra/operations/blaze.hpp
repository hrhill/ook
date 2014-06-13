#ifndef LINEAR_ALGEBRA_OPERATIONS_BLAZE_HPP_
#define LINEAR_ALGEBRA_OPERATIONS_BLAZE_HPP_

#if HAVE_BLAZE

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

} // ns linalg

#endif

#endif
