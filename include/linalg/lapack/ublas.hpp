#ifndef LINALG_LAPACK_UBLAS_HPP_
#define LINALG_LAPACK_UBLAS_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include <boost/numeric/bindings/blas/level1.hpp>
#include <boost/numeric/bindings/blas/level2.hpp>
#include <boost/numeric/bindings/blas/level3.hpp>

#include <boost/numeric/bindings/lapack/computational/potrf.hpp>
#include <boost/numeric/bindings/lapack/computational/potrs.hpp>
#include <boost/numeric/bindings/lapack/computational/potri.hpp>
#include <boost/numeric/bindings/lapack/computational/geqrf.hpp>
#include <boost/numeric/bindings/lapack/computational/tptri.hpp>
#include <boost/numeric/bindings/lapack/driver/posv.hpp>

#include <boost/numeric/bindings/std/vector.hpp>

#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>

#include <boost/numeric/bindings/trans.hpp>
#include <boost/numeric/bindings/symm.hpp>
#include <boost/numeric/bindings/column.hpp>

namespace linalg{

template <typename Matrix>
int
potrf(Matrix& a)
{
    boost::numeric::ublas::symmetric_adaptor<
        Matrix,
        boost::numeric::ublas::lower> sqrt_a(a);
    return boost::numeric::bindings::lapack::potrf(sqrt_a);
}

template <typename MatrixA, typename MatrixB>
int
potrs(MatrixA& a, MatrixB& b)
{
    boost::numeric::ublas::symmetric_adaptor<
        MatrixA,
        boost::numeric::ublas::lower> sa(a);
    return boost::numeric::bindings::lapack::potrs(sa, b);
}

template <typename MatrixA, typename MatrixB>
int
posv(MatrixA& a, MatrixB& b)
{
    boost::numeric::ublas::symmetric_adaptor<
        MatrixA,
        boost::numeric::ublas::lower> sa(a);
    return boost::numeric::bindings::lapack::posv(sa, b);
}

template <typename Matrix>
int
potri(Matrix& a)
{
    boost::numeric::ublas::symmetric_adaptor<
        Matrix,
        boost::numeric::ublas::lower> ia(a);
    return boost::numeric::bindings::lapack::potri(ia);
}

template <typename Matrix>
int
geqrf(Matrix& A, std::vector<int>& tau)
{
    return boost::numeric::bindings::lapack::geqrf(A, tau);
}

} // ns linalg

#endif
