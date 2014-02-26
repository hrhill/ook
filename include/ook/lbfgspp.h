#ifndef OOK_LINE_SEARCH_METHODS_LBFGSPP_H_
#define OOK_LINE_SEARCH_METHODS_LBFGSPP_H_

#include <tuple>
#include <boost/numeric/ublas/matrix.hpp>

#include "ook/line_search_method.h"

namespace ook{

namespace detail{

/// \brief Implementation of the required steps of line_search_method
/// for BFGS method.
template <typename X>
struct lbfgspp{
    typedef X vector_type;
    typedef typename X::value_type value_type;
    typedef state<X> state_type;

    template <typename F>
    static
    state_type
    initialise(F objective_function, const X& x0)
    {
        state_type s(x0.size(), true);
        std::tie(s.fx, s.dfx) = objective_function(x0);
        s.dfx0 = s.dfx;
        return s;
    }

    static
    vector_type
    descent_direction(state_type& s)
    {
        s.p = -boost::numeric::ublas::prod(s.H, s.dfx);
        return s.p;
    }

    static
    state_type
    update(state_type s){
        namespace ublas = boost::numeric::ublas;
        typedef ublas::matrix<double, ublas::column_major> matrix_type;

        X dx = s.a * s.p;
        X y(s.dfx - s.dfx0);
        const int n = s.dfx.size();

        const value_type rho = 1.0/detail::inner_product(y, dx);
        matrix_type Z(ublas::identity_matrix<double>(n) - rho * ublas::outer_prod(dx, y));
        matrix_type ss = rho * ublas::outer_prod(dx, dx);

        if (s.iteration == 1){
            const value_type hii = detail::inner_product(dx, dx);
            for (int i = 0; i < n; ++i){
                s.H(i, i) = hii;
            }
        }
        s.H = ublas::prod(Z, matrix_type(ublas::prod(s.H, ublas::trans(Z)))) + ss;
        s.dfx0 = s.dfx;
        return s;
    }
};

} // ns detail

/// \brief The Broyden-Fletcher-Goldfarb-Shanno (BFGS) algorithm.
/** \details Implementation of the BFGS algorithm using the generic line search function.
**/
template <typename F, typename X, typename Options, typename Stream>
std::tuple<ook::message, X>
bfgs(F objective_function, const X& x0, const Options& opts, Stream& stream)
{
    return line_search_method<detail::bfgs<X>>(objective_function, x0, opts, stream);
}

} //ns ook


#endif
