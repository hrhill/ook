# ifndef OOK_FACTORISATIONS_CHOLESKY_H_
# define OOK_FACTORISATIONS_CHOLESKY_H_

#include <vector>
#include <stdexcept>

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

#include <boost/numeric/bindings/std/vector.hpp>

#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>

namespace ook{

// select the lower triangular part of the matrix.
template <typename Matrix>
Matrix
select_lower_triangular(const Matrix& m){

    const int nrows = m.size1();
    const int ncols = m.size2();

    const int dim = std::min(nrows, ncols);
    const int row_offset = std::max(0, nrows - ncols);
    Matrix lower(dim, dim, 0.0);

    for (int i = 0; i < dim; ++i){
        for (int j = 0; j <= i; ++j){
            lower(i, j) = m(row_offset + i, j);
        }
    }
    return lower;
}

// select the upper triangular part of the matrix.
template <typename Matrix>
Matrix
select_upper_triangular(const Matrix& m){

    const int nrows = m.size1();
    const int ncols = m.size2();

    const int dim = std::min(nrows, ncols);
    const int col_offset = std::max(0, ncols - nrows);    
    Matrix upper(dim, dim, 0.0);

    for (int i = 0; i < dim; ++i){
        for (int j = i; j < dim; ++j){
            upper(i, j) = m(i, col_offset + j);
        }
    }
    return upper;
}

// get the lower triangular cholesky factor of m.
template <typename Matrix>
Matrix
get_lower_cholesky_factor(Matrix m){
    namespace ublas = boost::numeric::ublas;
    namespace lapack = boost::numeric::bindings::lapack;

    ublas::symmetric_adaptor<Matrix, ublas::lower> sqrt_m(m);
    int info = lapack::potrf(sqrt_m);
    if (info == -1){
        std::cout << "\n\nException: Trouble factoring " << sqrt_m << std::endl;
        throw std::runtime_error("get_lower_cholesky_factor lapack::potrf " + std::to_string(info));
    }

    return select_lower_triangular<Matrix>(sqrt_m);
}

template <typename MatrixType, typename VectorType>
VectorType
cholesky_solve(MatrixType A, const VectorType& y){

	namespace lapack = boost::numeric::bindings::lapack;
	namespace ublas = boost::numeric::ublas;
    /// Need to adapt types to that they are compatible with lapack
    ublas::symmetric_adaptor<MatrixType, ublas::lower> sa(A);

    // factorize
    int info = lapack::potrf(sa);
    if (info == -1)
        throw std::runtime_error("potrf failed in cholesky_solve" + std::to_string(info));

    // set up Ax = y
    MatrixType y1(y.size(), 1);
    ublas::column(y1, 0) = y;

    // solve
    lapack::potrs(sa, y1);

    return ublas::column(y1, 0);
}

template <typename MatrixType>
MatrixType
cholesky_invert(MatrixType A){

	namespace lapack = boost::numeric::bindings::lapack;
	namespace ublas = boost::numeric::ublas;

    /// Need to adapt types to that they are compatible with lapack
    ublas::symmetric_adaptor<MatrixType, ublas::lower> sa(A);

    // factorize
    int info = lapack::potrf(sa);
    if (info == -1){
        throw std::runtime_error("potrf failed in cholesky_invert " + std::to_string(info));
    }

    lapack::potri(sa);

    return sa;
}

template <typename MatrixType>
typename MatrixType::value_type
log_cholesky_determinant(MatrixType A){

	namespace lapack = boost::numeric::bindings::lapack;
	namespace ublas = boost::numeric::ublas;

    ublas::symmetric_adaptor<MatrixType, ublas::lower> sa(A);

    // factorize
    int info = lapack::potrf(sa);
    if (info == -1)
        throw std::runtime_error("potrf failed in cholesky_determinant" + std::to_string(info));

    double logd(0.0);
    for (typename MatrixType::size_type i = 0; i < A.size1(); ++i){
        logd += log(sa(i, i));
    }
    // return the square since |A| = |L|^2
    return  2.0 * logd;
}

template <typename MatrixType>
typename MatrixType::value_type
cholesky_determinant(MatrixType A){

    return  exp(log_cholesky_determinant(A));
}

template <typename MatrixType>
MatrixType
ldlt_factorisation(MatrixType A)
{
    typedef typename MatrixType::size_type size_type;
    const size_type n = A.size1();

    MatrixType L(n, n, 0.0);    
    MatrixType c(n, n, 0.0);
    for (size_type j = 0; j < n; ++j){
        c(j, j) = A(j, j);
        for (size_type s = 0; s < j; ++s){
            c(j, j) -= L(s, s) * std::pow(L(j, s), 2);
        }
        L(j, j) = c(j, j);

        for (size_type i = j + 1; i < n; ++i){
            c(i, j) = A(i, j);
            for (size_type s = 0; s < j; ++s){
                c(i, j) -= L(s, s) * L(i, s) * L(j, s);
            }
            L(i, j) = c(i, j) / L(j, j);
        }
    }
/*
    MatrixType L = get_lower_cholesky_factor(A);
    MatrixType D(n, n, 0.0);
    for (size_type i = 0; i < n; ++i){
        D(i, i) = sqrt(L(i, i));
        L(i, i) = 1.0;
    }
*/    
    return L;
    //return boost::numeric::ublas::prod(L, D);
}

} // ns ook

#endif
