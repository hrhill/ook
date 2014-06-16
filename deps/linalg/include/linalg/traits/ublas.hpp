#ifndef LINALG_TRAITS_UBLAS_HPP_
#define LINALG_TRAITS_UBLAS_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace linalg{

template <typename T>
struct associated_matrix;

template <typename T>
struct associated_matrix<boost::numeric::ublas::vector<T>>
{
    using type = boost::numeric::ublas::matrix<T, boost::numeric::ublas::column_major>;
};


} // ns linalg

#endif
