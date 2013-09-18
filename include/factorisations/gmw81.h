#ifndef MODIFIED_LDLT_H_
#define MODIFIED_LDLT_H_

#include <cassert>
#include <algorithm>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

namespace ook{

namespace detail{
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


template <typename Vector>
std::tuple<typename Vector::size_type, typename Vector::value_type>
max_magnitude(const Vector& x)
{
    typedef typename Vector::size_type size_type;
    typedef typename Vector::value_type value_type;  

    const size_type n = x.size();

    if (n == 0)
        return std::make_tuple(0, 0);
    size_type idx = 0;
    value_type xmax = fabs(x(0));
    for (size_type i = 0; i < n; ++i){
        const value_type xi = fabs(x(i));
        idx = xmax < xi ?  i : idx;
        xmax = idx == i ? xi : xmax;
    }
    return std::make_tuple(idx, xmax);      
}

} // ns detail

/// \brief Implementation of the modified cholesky factorisation
/// of Gill, Murray and Wright from Ch.4 page 111 of "Practical Optimisation".
/// The algorithm is commonly referred to as gmw81, hence the name here.
/// The function takes a symmetric matrix G, and returns a modified variant
/// of the LDLt factorisation. The returned matrix has layout,
/// \f[
/// \begin{pmatrix}
///  d_{11} & 0      & \cdots  & 0 \\
///  l_{21} & d_{22} & \cdots  & \vdots\\
///  \vdots &        & \ddots  & \vdots\\
///  l_{nn} & l_{n2} & \cdots   & d_{nn}\\
/// \end{pmatrix}
/// \f] 
///  
template <typename Matrix>
Matrix
gmw81(Matrix G)
{
    using namespace boost::numeric::ublas;

    typedef typename Matrix::value_type real_type;
    typedef typename Matrix::size_type size_type;    

    // MC1
    const size_type n = G.size1();
    const real_type eps = std::numeric_limits<real_type>::epsilon();
    const real_type nu = std::max(1.0, sqrt(std::pow(n, 2) - 1.0));
    const real_type gamma = std::get<1>(detail::max_magnitude_diagonal(G));
    const real_type eta = detail::max_magnitude_off_diagonal(G);    
    const real_type beta2 = std::max({gamma, eta/nu, eps});
    const real_type delta = sqrt(eps);

    // MC2
    Matrix c(n, n, 0);
    Matrix L(n, n, 0);
    Matrix D(n, 1, 0);
    Matrix E(n, 1, 0);

    for (size_type i = 0; i < n; ++i){
        c(i, i) = G(i, i);
    }
    for (int j = 0; j < n; ++j){
        // MC3 
        real_type cqq;
        size_type q;
        std::tie(q, cqq) = detail::max_magnitude_diagonal(matrix_range<Matrix>(c, range(j, n), range(j, n)));

        // Pivoting
        if (q != j){
            matrix_row<Matrix> rowq(G, q);
            matrix_row<Matrix> rowj(G, j);
            rowq.swap(rowj);

            matrix_row<Matrix> colq(G, q);
            matrix_row<Matrix> colj(G, j);
            colq.swap(colj);
        }
        // MC4
        for (int s = 0; s < j; ++s){
            L(j, s) = c(j, s)/L(s, 0);
        }
        for (int i = j + 1; i < n; ++i){
            c(i, j) = G(i, j);
            for (int s = 0; s < j; ++s){
                c(i, j) -= L(j, s) * c(i, s);
            }
        }
        size_type theta_idx;
        real_type theta = 0;
        if (j < n - 1)
            std::tie(theta_idx, theta) = detail::max_magnitude(matrix_column<Matrix>(c, j));

        // MC 5
        L(j, j) = std::max({delta, fabs(c(j, j)), std::pow(theta, 2)/beta2});

        for (int i = j + 1; i < n; ++i){
            c(i, i) -= std::pow(c(i, j), 2)/L(j, j);
        }
        E(j, 0) = L(j, j) - c(j, j);
    }
    return L;
}

} // ns ook

#endif
