#ifndef LINEAR_ALGEBRA_OPERATIONS_UBLAS_HPP_
#define LINEAR_ALGEBRA_OPERATIONS_UBLAS_HPP_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace linalg{

template <
    template<typename, typename> class V,
    typename T,
    typename A>
size_t
size(const V<T, A>& v)
{
    return v.size();
}

template <
    template <typename, typename, typename> class Matrix,
    typename T,
    typename F,
    typename A
>
size_t
num_rows(const Matrix<T, F, A>& m)
{
    return m.size1();
}

template <typename M>
size_t
num_rows(const boost::numeric::ublas::matrix_range<M>& m){
    return m.size1();
}

template <
    template <typename, typename, typename> class Matrix,
    typename T,
    typename F,
    typename A
>
size_t
num_cols(const Matrix<T, F, A>& m)
{
    return m.size2();
}

template <typename M>
size_t
num_cols(const boost::numeric::ublas::matrix_range<M>& m){
    return m.size1();
}


template <typename Matrix>
auto
column(Matrix& a, size_t idx) -> decltype(boost::numeric::ublas::column(a, idx))
{
    return boost::numeric::ublas::column(a, idx);
}

template <typename Matrix>
auto
row(Matrix& a, size_t idx) -> decltype(boost::numeric::ublas::row(a, idx))
{
    return boost::numeric::ublas::row(a, idx);
}

template <typename MatrixType>
inline
MatrixType
trans(const MatrixType& M)
{
	return boost::numeric::ublas::trans(M);
}

template <typename Vector>
auto
outer_prod(const Vector& x, const Vector& y) -> decltype(boost::numeric::ublas::outer_prod(x, y))
{
    return boost::numeric::ublas::outer_prod(x, y);
}

template <typename T, typename U>
void swap(T& t, U& u){
    std::swap(t, u);
}

} // ns linalg

#endif
