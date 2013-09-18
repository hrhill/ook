#ifndef MODIFIED_LDLT_H_
#define MODIFIED_LDLT_H_

#include <cassert>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

template <class Matrix>
Matrix
cholesky(const Matrix& A)
{
    typedef typename Matrix::value_type T;

    const size_t n = A.size1();
    Matrix L(n, n, 0);

    for (size_t j = 0; j < n; ++j){
        // D entry
        T cjj = A(j, j);
        for (size_t s = 0; s < j; ++s){
            cjj -= std::pow(L(j, s), 2);
        }
        L(j, j) = sqrt(cjj);

        // Lower diagonal
        for (size_t i = j + 1; i < n; ++i){
            T cij = A(i, j);
            for (size_t s = 0; s < j; ++s){
                cij -= L(i, s) * L(j, s);
            }
            L(i, j) = cij/L(j, j);
        }
    }
    return L;
}

template <class Matrix>
Matrix
cholesky_ldlt(const Matrix& A)
{
    typedef typename Matrix::value_type T;

    const size_t n = A.size1();
    Matrix L(n, n, 0);

    for (size_t j = 0; j < n; ++j){
        // D entry
        T cjj = A(j, j);
        for (size_t s = 0; s < j; ++s){
            cjj -= L(s, s) * std::pow(L(j, s), 2);
        }
        L(j, j) = cjj;

        // Lower diagonal
        for (size_t i = j + 1; i < n; ++i){
            T cij = A(i, j);
            for (size_t s = 0; s < j; ++s){
                cij -= L(s, s) * L(i, s) * L(j, s);
            }
            L(i, j) = cij / L(j, j);
        }
    }
    return L;
}

#endif
