#ifndef OOK_LINE_SEARCH_METHODS_NEWTON_H_
#define OOK_LINE_SEARCH_METHODS_NEWTON_H_

#include <tuple>

#include "line_search_method.h"

#include "factorisations/gmw81.h"

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

template <typename F, typename X, typename Options>
std::tuple<ook::state_value, X>
newton(F objective_function, const X& x0, const Options& opts)
{
    return line_search_method<detail::newton<X>>(objective_function, x0, opts);
}

} //ns ook


#endif
