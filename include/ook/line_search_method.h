#ifndef OOK_LINE_SEARCH_METHOD_H_
#define OOK_LINE_SEARCH_METHOD_H_

#include <iostream>
#include <iomanip>
#include <boost/numeric/ublas/matrix.hpp>

#include "ook/norms.h"
#include "ook/message.h"

#include "ook/line_search/more_thuente/more_thuente.h"

namespace ook{

namespace detail{

enum class state_tag{init, iterate, final};

template <typename X>
struct
state{
    typedef X vector_type;
    typedef typename X::value_type value_type;
    typedef boost::numeric::ublas::matrix<value_type,
              boost::numeric::ublas::column_major> matrix_type;

    state(const int n, const bool with_matrix = false)
    :
         dfx(n),
         dfx0(n),
         p(n),
         dx(n),
         H(n * with_matrix, n * with_matrix, 0.0),
         a(1),
         beta(0),
         iteration(0),
         tag(state_tag::init)
    {
        if (with_matrix){
            for (int i = 0; i < n; ++i){
                H(i, i) = 1.0;
            }
        }
    }

    friend
    std::ostream&
    operator<<(std::ostream& out, const state& s){
        if (s.tag == state_tag::init){
            out << std::endl
                << std::setw(6) << "n"
                << std::setw(6) << "nfev"
                << std::scientific
                << std::setw(14) << "a"
                << std::setw(14) << "fx"
                << std::setw(14) << "max ||dfx||"
                << std::setw(14) << "max ||dx||" << std::endl;
        }

        if (s.tag == state_tag::iterate){
            out << std::setw(6) << s.iteration
                << std::setw(6) << 0//s.nfev_total
                << std::scientific
                << std::setw(14) << s.a
                << std::setw(14) << s.fx
                << std::setw(14) << ook::norm_infinity(s.dfx)
                << std::setw(14) << ook::norm_infinity(s.dx);
        }

        if (s.tag == state_tag::final)
        {
            out << "\nstatus : " << s.msg << std::endl;
            out << std::setw(8) << "iter"
                << std::setw(8) << "nfev"
                << std::setw(16) << "fx"
                << std::setw(16) << "max ||dfx||"
                << std::setw(16) << "max ||dx||" << std::endl;

            out << std::setw(8) << s.iteration
                << std::setw(8) << 0//s.nfev_total
                << std::scientific
                << std::setw(16) << s.fx
                << std::setw(16) << ook::norm_infinity(s.dfx)
                << std::setw(16) << ook::norm_infinity(s.dx) << std::endl;
        }
        return out;
    }

    value_type fx;
    vector_type dfx;
    vector_type dfx0;
    vector_type p;
    vector_type dx;
    matrix_type H;
    value_type a;
    value_type beta;
    int iteration;
    state_tag tag;
    message msg;
};

template <typename T>
typename T::value_type
inner_product(const T& x, const T& y){
    typedef typename T::value_type value_type;
    return std::inner_product(x.begin(), x.end(), y.begin(), value_type(0.0));
}

// Meta program to select the right function call
// based on the properties of the return type.
template <typename F, typename X, typename State, int dim>
struct function_caller{};

template <typename F, typename X, typename State>
struct function_caller<F, X, State, 2>{
    static
    void
    call(F objective_function, const X& x, State& s){
        std::tie(s.fx, s.dfx) = objective_function(x);
    }
};

template <typename F, typename X, typename State>
struct function_caller<F, X, State, 3>{
    static
    void
    call(F objective_function, const X& x, State& s){
        std::tie(s.fx, s.dfx, s.H) =  objective_function(x);
    }
};


} // ns detail

template <typename Scheme, typename F, typename X, typename Options, typename Observer>
std::tuple<ook::message, X>
line_search_method(F objective_function, X x, const Options& opts, Observer& observer)
{
    typedef typename X::value_type real_type;
    typedef typename Scheme::state_type state_type;

    typedef decltype(objective_function(x)) result_type;
    typedef detail::function_caller<F, X, state_type, std::tuple_size<result_type>::value> fcaller_type;

    const real_type epsilon = std::numeric_limits<real_type>::epsilon();
    const real_type dx_eps = sqrt(epsilon);
    const real_type df_eps = exp(log(epsilon)/3);

    state_type s = Scheme::initialise(objective_function, x);
    X dx(x.size());
    s.iteration = 0;
    uint nfev_total = 0;

    observer(s);
    do {
        s.tag = detail::state_tag::iterate;
        // Get descent direction and set up line search procedure.
        X p = Scheme::descent_direction(s);
        real_type dfx_dot_p = detail::inner_product(s.dfx, p);
        // do line search
        uint nfev = 0;
        s.a = 1.0;
        // take a reference to the state variable, ensuring that fx and dfx get updated
        auto phi = [&nfev, &s, &dfx_dot_p, &x, &p, objective_function](const real_type& a){
            ++nfev;
            fcaller_type::call(objective_function, x + a * p, s);
            dfx_dot_p = detail::inner_product(s.dfx, p);
            return std::make_pair(s.fx, dfx_dot_p);
        };

        // Store current fx value since line search overwrites the state values.
        const real_type fxk = s.fx;
        std::tie(s.msg, s.a) = ook::line_search::more_thuente(phi, s.fx, dfx_dot_p, s.a, opts);

        if (s.msg != ook::message::convergence){
            break;
        }
        dx = s.a * p;
        x += dx;
        nfev_total += nfev;
        ++s.iteration;

        // Convergence criteria assessment base on p306 in Gill, Murray and Wright.
        const double theta = epsilon * (1 + fabs(s.fx));
        const bool u1 = (fxk - s.fx) <= theta;
        const bool u2 = ook::norm_infinity(dx) <=  dx_eps * (1.0 + ook::norm_infinity(x));
        const bool u3 = ook::norm_infinity(s.dfx) <= df_eps * (1.0 + fabs(s.fx));

        observer(s);
        if ((u1 & u2) || u3){
            s.msg = ook::message::convergence;
            break;
        }

        s = Scheme::update(s);

    } while(true);

    s.tag = detail::state_tag::final;
    observer(s);
    return std::make_pair(s.msg, x);
}

} // ns ook
#endif
