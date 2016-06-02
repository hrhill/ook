#ifndef OOK_OOK_MATRIX_HPP_
#define OOK_OOK_MATRIX_HPP_

#include <blaze/Math.h>

namespace ook
{
inline namespace v1
{

typedef blaze::DynamicMatrix<double> matrix;
typedef blaze::CompressedMatrix<double> sparse_matrix;
typedef blaze::SymmetricMatrix<matrix> symmetric_matrix;
typedef blaze::SymmetricMatrix<sparse_matrix> sparse_symmetric_matrix;

matrix eye(size_t n);
matrix sympd(size_t n, std::mt19937&);
matrix n01_matrix(size_t m, size_t n, std::mt19937&);

double norm_inf(const matrix& a);
double norm_inf(const symmetric_matrix& a);

}
}

#endif