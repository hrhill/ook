# ifndef OOK_FACTORISATIONS_CHOLESKY_H_
# define OOK_FACTORISATIONS_CHOLESKY_H_

#include <vector>
#include <stdexcept>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/triangular.hpp>

#include <boost/numeric/bindings/blas/level1.hpp>
#include <boost/numeric/bindings/blas/level2.hpp>
#include <boost/numeric/bindings/blas/level3.hpp>

#include <boost/numeric/bindings/lapack/computational/potrf.hpp>
#include <boost/numeric/bindings/lapack/computational/potrs.hpp>
#include <boost/numeric/bindings/lapack/computational/potri.hpp>
#include <boost/numeric/bindings/lapack/computational/geqrf.hpp>
#include <boost/numeric/bindings/lapack/computational/tptri.hpp>

#include <boost/numeric/bindings/std/vector.hpp>

#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>

#include "./tools.h"

namespace ook{
namespace factorisations{

// get the lower triangular cholesky factor of m.
template <typename Matrix>
Matrix
cholesky(Matrix m){
    namespace ublas = boost::numeric::ublas;
    namespace lapack = boost::numeric::bindings::lapack;

    ublas::symmetric_adaptor<Matrix, ublas::lower> sqrt_m(m);
    int info = lapack::potrf(sqrt_m);
    if (info == -1){
        std::cout << "\n\nException: Trouble factoring " << sqrt_m << std::endl;
        throw std::runtime_error("get_lower_cholesky_factor lapack::potrf " + std::to_string(info));
    }

    return tools::select_lower_triangular<Matrix>(sqrt_m);
}

}
} // ns ook

#endif
