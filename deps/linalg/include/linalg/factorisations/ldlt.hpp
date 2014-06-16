# ifndef LINALG_FACTORISATIONS_LDLT_HPP_
# define LINALG_FACTORISATIONS_LDLT_HPP_

#include <vector>
#include <stdexcept>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include "linalg/factorisations/tools.hpp"

namespace linalg{
namespace factorisations{

/// \brief LDL^t factorisation.
template <typename Matrix>
Matrix
ldlt(Matrix A)
{
    using namespace boost::numeric::ublas;

    const int n = linalg::num_rows(A);

    Matrix L(n, n, 0.0);
    Matrix c(n, n, 0.0);
    for (size_t j = 0; j < n; ++j){
/*
        real_type cqq;
        size_t q;
        std::tie(q, cqq) = tools::max_magnitude_diagonal(matrix_range<Matrix>(c, range(j, n), range(j, n)));

        // Pivoting
        if (q != j){
            matrix_row<Matrix> rowq(A, q);
            matrix_row<Matrix> rowj(A, j);
            rowq.swap(rowj);

            matrix_column<Matrix> colq(A, q);
            matrix_column<Matrix> colj(A, j);

            colq.swap(colj);
        }
*/
        c(j, j) = A(j, j);
        for (size_t s = 0; s < j; ++s){
            c(j, j) -= L(s, s) * std::pow(L(j, s), 2);
        }
        L(j, j) = c(j, j);

        for (size_t i = j + 1; i < n; ++i){
            c(i, j) = A(i, j);
            for (size_t s = 0; s < j; ++s){
                c(i, j) -= L(s, s) * L(i, s) * L(j, s);
            }
            L(i, j) = c(i, j) / L(j, j);
        }
    }
    return L;
}

} // ns factorisations

} // ns ook

#endif
