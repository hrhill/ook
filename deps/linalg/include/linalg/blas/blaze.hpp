#ifndef LINALG_BLAS_BLAZE_HPP_
#define LINALG_BLAS_BLAZE_HPP_

#if HAVE_BLAZE

#include <blaze/Math.h>

namespace linalg{

template <
    template <typename, bool> class MatrixType,
    typename T1,
    bool SO,
    template <typename, bool> class VectorType,
    typename T2,
    bool TF>
inline
void
gemv(const T1& a,
     const MatrixType<T1, SO>& m,
     const VectorType<T2, TF>& v,
     const T1& b,
     VectorType<T2, TF>& res)
{
    res = b * res;
    res += a * m * v;
}

/// gemm
template <
    template <typename, bool> class MatrixType,
    typename T,
    bool SO>
inline
void
gemm(const T& a,
        const MatrixType<T, SO>& A,
        const MatrixType<T, SO>& B,
        const T& b,
        MatrixType<T, SO>& C)
{
    C = b * C;
    C += a * A * B;
}

} //ns linalg

#endif

#endif
