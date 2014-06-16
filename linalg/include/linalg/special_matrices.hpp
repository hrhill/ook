#ifndef LINALG_SPECIAL_MATRICES_HPP_
#define LINALG_SPECIAL_MATRICES_HPP_

namespace linalg{

template <typename MatrixType>
MatrixType
constant_diagonal_matrix(const int n, const double d)
{
    MatrixType m(n, n, 0.0);
    for (int i = 0; i < n; ++i){
        m(i, i) = d;
    }
    return m;
}

template <typename MatrixType>
MatrixType
identity_matrix(const int n)
{
    return constant_diagonal_matrix<MatrixType>(n, 1.0);
}

} // ns linalg

#endif
