#ifndef LINEAR_ALGEBRA_TRAITS_BLAZE_HPP_
#define LINEAR_ALGEBRA_TRAITS_BLAZE_HPP_

#if HAVE_BLAZE

#include <blaze/Math.h>

namespace linalg{

template <T>
struct associated_matrix_type;

template <T>
struct associated_matrix_type<blaze::DynamicVector<double>>
{
    using matrix_type = blaze::DynamicMatrix<double, blaze::columnMajor>;
};

} // ns linalg

#endif

#endif
