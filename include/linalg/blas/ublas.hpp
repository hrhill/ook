#ifndef LINALG_BLAS_UBLAS_HPP_
#define LINALG_BLAS_UBLAS_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include <boost/numeric/bindings/blas/level1.hpp>
#include <boost/numeric/bindings/blas/level2.hpp>
#include <boost/numeric/bindings/blas/level3.hpp>

#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>

namespace linalg{

template <typename Matrix, typename Vector>
inline
void
gemv(const typename Vector::value_type& a,
        const Matrix& m,
        const Vector& v,
        const typename Vector::value_type& b,
        Vector& res)
{
    boost::numeric::bindings::blas::gemv(a, m, v, b, res);
}

template <typename Matrix1, typename Matrix2, typename Matrix3>
inline
void
gemm(const typename Matrix1::value_type& a,
     const Matrix1& A,
     const Matrix2& B,
     const typename Matrix1::value_type& b,
     Matrix3& C)
{
    boost::numeric::bindings::blas::gemm(a, A, B, b, C);
}



} // ns linalg

#endif
