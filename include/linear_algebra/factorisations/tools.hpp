#ifndef LINEAR_ALGEBRA_FACTORISATIONS_TOOLS_H_
#define LINEAR_ALGEBRA_FACTORISATIONS_TOOLS_H_

#include <algorithm>

namespace linalg{
namespace factorisations{
namespace tools{

///\brief Select the lower triangular part of the matrix.
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

/// \brief Select the upper triangular part of the matrix.
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

/// \brief Get the value of maximum magnitude on the diagonal of the matrix.
template <typename Matrix>
std::tuple<typename Matrix::size_type, typename Matrix::value_type>
max_magnitude_diagonal(const Matrix& m){

    typedef typename Matrix::size_type size_type;
    typedef typename Matrix::size_type value_type;

    const size_type n = m.size1();

    if (n == 0)
        return std::make_tuple(0, 0);

    size_type idx = 0;
    value_type mmd = fabs(m(0, 0));
    for (size_type i = 1; i < n; ++i){
        const value_type mii = fabs(m(i, i));
        idx = mmd < mii ?  i : idx;
        mmd = idx == i ? mii : mmd;
    }
    return std::make_tuple(idx, mmd);
}

/// \brief Get the value of maximum magnitude off the diagonal of the matrix.
template <typename Matrix>
typename Matrix::value_type
max_magnitude_off_diagonal(const Matrix& m){

    typedef typename Matrix::size_type size_type;
    typedef typename Matrix::value_type value_type;

    const size_type n = m.size1();

    value_type mmd = fabs(m(0, 0));
    for (size_type i = 0; i < n; ++i){
        for (size_type j = i + 1; j < n; ++j){
            mmd = std::max(fabs(m(i, j)), mmd);
        }
    }
    return mmd;
}

} // ns tools
} // ns factorisations
} // ns ook

#endif
