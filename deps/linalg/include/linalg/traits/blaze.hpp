#ifndef LINALG_TRAITS_BLAZE_HPP_
#define LINALG_TRAITS_BLAZE_HPP_

#if HAVE_BLAZE

#include <blaze/Math.h>

namespace linalg{

template <typename T>
struct associated_matrix;

template <typename T>
struct associated_matrix<blaze::DynamicVector<T>>
{
    using type = blaze::DynamicMatrix<T, blaze::columnMajor>;
};

} // ns linalg

#endif

#endif
