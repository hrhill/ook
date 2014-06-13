#ifndef LINEAR_ALGEBRA_TRAITS_UBLAS_HPP_
#define LINEAR_ALGEBRA_TRAITS_UBLAS_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

namespace linalg{

template <T>
struct associated_matrix_type;

template <T>
struct associated_matrix_type<boost::numeric::ublas::vector<T>>
{
    using matrix_type = boost::numeric::ublas::matrix<T, boost::numeric::column_major>;
};


} // ns linalg

#endif
