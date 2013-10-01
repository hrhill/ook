#ifndef OOK_LINE_SEARCH_METHODS_NEWTON_H_
#define OOK_LINE_SEARCH_METHODS_NEWTON_H_

#include <tuple>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/assignment.hpp>

#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/lapack/computational/potrs.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/symmetric.hpp>

#include "ook/line_search_methods/line_search_method.h"
#include "ook/factorisations/gmw81.h"

namespace ook{

namespace detail{

template <typename Matrix>
Matrix
convert_to_cholesky(const Matrix& LD)
{
    int n = LD.size1();
    Matrix L(LD);
    Matrix D(n, n, 0);

    for (int i = 0; i < n; ++i){
        D(i, i) = sqrt(L(i, i));
        L(i, i) = 1.0;
    }
    return boost::numeric::ublas::prod(L, D);
}

template <typename Matrix, typename Vector>
Vector
solve(Matrix A, const Vector& b)
{
    Matrix LD = ook::factorisations::gmw81(A);
    Matrix L = convert_to_cholesky(LD);

    boost::numeric::ublas::symmetric_adaptor<Matrix, boost::numeric::ublas::lower> sa(L);    
    Matrix b1(b.size(), 1);

    boost::numeric::ublas::column(b1, 0) = b;
    boost::numeric::bindings::lapack::potrs(sa, b1);

    return boost::numeric::ublas::column(b1, 0);
}

/// \brief Implementation of the required steps of line_search_method
/// for Newtons method.
template <typename X>
struct newton{
    typedef X vector_type;
    typedef typename X::value_type value_type;
    typedef state<X> state_type;

    template <typename F>
    static
    state_type
    initialise(F objective_function, const X& x0)
    {
        state_type s(x0.size());
        std::tie(s.fx, s.dfx, s.d2fx) = objective_function(x0);
        return s;
    }
    
    static
    vector_type
    descent_direction(state_type& s)
    {
        return -detail::solve(s.d2fx, s.dfx);
    }

    static
    state_type
    update(state_type s){
        return s;
    }    
};

} // ns detail

/// \brief The Newton algorithm.
/** \details Implementation of the Newton algorithm using the generic line search function.
**/
template <typename F, typename X, typename Options, typename Stream>
std::tuple<ook::state_value, X>
newton(F objective_function, const X& x0, const Options& opts, Stream& stream)
{
    return line_search_method<detail::newton<X>>(objective_function, x0, opts, stream);
}

} //ns ook


#endif
