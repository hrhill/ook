#ifndef LINALG_FACTORISATIONS_GMW81_HPP_
#define LINALG_FACTORISATIONS_GMW81_HPP_

#include <cassert>
#include <algorithm>

#include "linalg/std_traits.hpp"
#include "linalg/operations.hpp"
#include "linalg/factorisations/tools.hpp"

namespace linalg{
namespace factorisations{
namespace detail{

template <typename Matrix>
auto
calculate_theta(const Matrix& c, const int j)
{
    const int n = linalg::num_rows(c);
    auto theta(0.0);
    if (j < n ){
        for (int i = j + 1; i < n; ++i){
            theta = std::max(theta, c(i, j));
        }
    }
    return theta;
}

} // ns detail

/// \brief Modified cholesky factorisation of Gill, Murray and Wright.
/**
\details Modified cholesky factorisation of Gill, Murray and Wright
from Ch.4 page 111 of "Practical Optimisation".
The algorithm is commonly referred to as gmw81, hence the name here.
The function takes a symmetric matrix G, and returns a modified variant
of the LDLt factorisation. The returned matrix has layout,
\f[
\begin{pmatrix}
    d_{11} & 0      & \cdots  & 0 \\
    l_{21} & d_{22} & \cdots  & \vdots\\
    \vdots &        & \ddots  & \vdots\\
    l_{nn} & l_{n2} & \cdots   & d_{nn}\\
    \end{pmatrix}
\f]
\tparam Matrix Templated by matrix type, but currently restricted to Boost.Ublas.
**/
template <typename Matrix>
Matrix
gmw81(Matrix G)
{
    using namespace boost::numeric::ublas;

    typedef  remove_const_reference<decltype(G(0,0))> value_type;

    // MC1
    const size_t n = linalg::num_rows(G);
    const value_type eps = std::numeric_limits<value_type>::epsilon();
    const value_type nu = std::max(1.0, sqrt(std::pow(n, 2) - 1.0));
    const value_type gamma = std::get<1>(tools::max_magnitude_diagonal(G));
    const value_type eta = tools::max_magnitude_off_diagonal(G);
    const value_type beta2 = std::max({gamma, eta/nu, eps});
    const value_type delta = sqrt(eps);

    // MC2
    Matrix c(n, n, 0);
    Matrix L(n, n, 0);

    for (size_t i = 0; i < n; ++i){
        c(i, i) = G(i, i);
    }
    for (size_t j = 0; j < n; ++j){
        // MC3
/*
        value_type cqq;
        size_t q;
        std::tie(q, cqq) = tools::max_magnitude_diagonal(matrix_range<Matrix>(c, range(j, n), range(j, n)));

        // Pivoting
        if (q != j){
            matrix_row<Matrix> rowq(G, q);
            matrix_row<Matrix> rowj(G, j);
            rowq.swap(rowj);

            matrix_column<Matrix> colq(G, q);
            matrix_column<Matrix> colj(G, j);
            colq.swap(colj);
        }
*/
        // MC4
        for (size_t s = 0; s < j; ++s){
            L(j, s) = c(j, s)/L(s, s);
        }
        for (size_t i = j + 1; i < n; ++i){
            c(i, j) = G(i, j);
            for (int s = 0; s < j; ++s){
                c(i, j) -= L(j, s) * c(i, s);
            }
        }
        value_type theta = detail::calculate_theta(c, j);
        // MC 5
        L(j, j) = std::max({delta, fabs(c(j, j)), std::pow(theta, 2)/beta2});

        for (size_t i = j + 1; i < n; ++i){
            c(i, i) -= std::pow(c(i, j), 2)/L(j, j);
        }
    }
    return L;
}

} // ns factorisation

} // ns ook

#endif
