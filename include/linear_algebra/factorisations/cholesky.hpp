# ifndef LINEAR_ALGEBRA_FACTORISATIONS_CHOLESKY_H_
# define LINEAR_ALGEBRA_FACTORISATIONS_CHOLESKY_H_

#include "linear_algebra/lapack.hpp"
#include "linear_algebra/factorisations/tools.hpp"

namespace linalg{
namespace factorisations{

/// \brief User friendly Cholesky factorisation.
template <typename Matrix>
Matrix
cholesky(Matrix m){
    linalg::potrf(m);
    return tools::select_lower_triangular<Matrix>(m);
}

}
} // ns ook

#endif
