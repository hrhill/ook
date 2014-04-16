// Copyright 2013 Harry Hill
//
// This file is part of ook.
//
// ook is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ook is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with ook.  If not, see <http://www.gnu.org/licenses/>.

# ifndef OOK_FACTORISATIONS_CHOLESKY_H_
# define OOK_FACTORISATIONS_CHOLESKY_H_

#include <vector>
#include <stdexcept>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <boost/numeric/bindings/lapack/computational/potrf.hpp>
#include <boost/numeric/bindings/lapack/computational/potrs.hpp>

#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>

#include "ook/factorisations/tools.h"

namespace ook{
namespace factorisations{

/// \brief User friendly Cholesky factorisation.
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
