#ifndef OOK_LINE_SEARCH_METHOD_H_
#define OOK_LINE_SEARCH_METHOD_H_

#include <iomanip>

#include <boost/numeric/ublas/matrix.hpp>

#include "ook/norms.h"
#include "ook/line_search_methods/message.h"

#include "ook/line_search_methods/more_thuente/more_thuente.h"

namespace ook{

namespace detail{

template <typename X>
struct 
state{
    typedef X vector_type;
    typedef typename X::value_type value_type;
    
    typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> matrix_type;
    state(const int n, const bool with_matrix = false)
    :
         dfx(n), 
         dfx0(n), 
         p(n), 
         a(1), 
         beta(0), 
         iteration(0)
    {
        if (with_matrix){
            H.resize(n, n);
            for (int i = 0; i < n; ++i){
                H(i, i) = 1.0;              
                for (int j = 0; j < i; ++j){
                    H(i, j) = 0.0;
                    H(j, i) = 0.0;                    
                }
            }
        }
    }

    value_type fx;
    vector_type dfx;
    vector_type dfx0;        
    vector_type p;    
    matrix_type H;
    value_type a;    
    value_type beta;
    int iteration; 
};

template <typename T>
typename T::value_type
inner_product(const T& x, const T& y){
    typedef typename T::value_type value_type;
    return std::inner_product(x.begin(), x.end(), y.begin(), value_type(0.0));
}

template <typename Stream, typename X>
void
report(Stream& stream, message msg, int iteration, int nfev_total, int nfev, double a, double fx, const X& dfx, const X& dx)
{
    stream << std::setw(6) << iteration
              << std::setw(6) << nfev_total
              << std::scientific 
              << std::setw(14) << a
              << std::setw(14) << fx
              << std::setw(14) << ook::norm_infinity(dfx)
              << std::setw(14) << ook::norm_infinity(dx) << std::endl;  
}

template <typename Stream, typename X>
void
final_report(Stream& stream, message msg, int iteration, int nfev_total,double fx, const X& dfx, const X& dx)
{
    stream << "status : " << msg << std::endl;
    stream << std::setw(8) << "iter"
           << std::setw(8) << "nfev"
           << std::setw(16) << "fx"
           << std::setw(16) << "max ||dfx||"
           << std::setw(16) << "max ||dx||" << std::endl;  
    stream << std::setw(8) << iteration 
           << std::setw(8) << nfev_total
           << std::scientific 
           << std::setw(16) << fx
           << std::setw(16) << ook::norm_infinity(dfx)
           << std::setw(16) << ook::norm_infinity(dx) << std::endl;  
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

template <typename Scheme, typename F, typename X, typename Options, typename Stream>
std::tuple<ook::message, X>
line_search_method(F objective_function, X x, const Options& opts, Stream& stream)
{
    typedef typename X::value_type real_type;
    typedef typename Scheme::state_type state_type;

    typedef decltype(objective_function(x)) result_type;
    typedef detail::function_caller<F, X, state_type, std::tuple_size<result_type>::value> fcaller_type;

    const real_type epsilon = std::numeric_limits<real_type>::epsilon();
    const real_type dx_eps = sqrt(epsilon);
    const real_type df_eps = exp(log(epsilon)/3);
    const real_type fx_max = std::numeric_limits<real_type>::max();
    const real_type fx_min = std::numeric_limits<real_type>::min();    

    state_type s = Scheme::initialise(objective_function, x);
    X dx(x.size());
    s.iteration = 0;
    uint nfev_total = 0;
    ook::message msg;

    stream << std::endl;
    stream << std::setw(6) << "n"
           << std::setw(6) << "nfev"
           << std::scientific 
           << std::setw(14) << "a"
           << std::setw(14) << "fx"
           << std::setw(14) << "max ||dfx||"
           << std::setw(14) << "max ||dx||" << std::endl;

    do {
        // Get descent direction and set up line search procedure.
        X p = Scheme::descent_direction(s);
        real_type dfx_dot_p = detail::inner_product(s.dfx, p);
        dfx_dot_p = std::min(dfx_dot_p, std::numeric_limits<real_type>::max());
        dfx_dot_p = std::max(dfx_dot_p, std::numeric_limits<real_type>::lowest());

        // do line search
        uint nfev = 0;        
        s.a = 1.0;
        // take a reference to the state variable, ensuring that fx and dfx get updated
        auto phi = [&nfev, &s, &dfx_dot_p, &x, &p, objective_function](const real_type& a){
            ++nfev;
            fcaller_type::call(objective_function, x + a * p, s);
            dfx_dot_p = detail::inner_product(s.dfx, p);
            s.fx = std::min(s.fx, std::numeric_limits<real_type>::max());
            s.fx = std::max(s.fx, std::numeric_limits<real_type>::lowest());
            dfx_dot_p = std::min(dfx_dot_p, std::numeric_limits<real_type>::max());
            dfx_dot_p = std::max(dfx_dot_p, std::numeric_limits<real_type>::lowest());
            return std::make_pair(s.fx, dfx_dot_p);
        };

        // Store current fx value since line search overwrites the state values.
        const real_type fxk = s.fx;
        std::tie(msg, s.a) = ook::line_search::more_thuente(phi, s.fx, dfx_dot_p, s.a, opts);

        if (msg != ook::message::convergence){
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

        detail::report(stream, msg, s.iteration, nfev_total, nfev, s.a, s.fx, s.dfx, dx);
        if ((u1 & u2) || u3){
            msg = ook::message::convergence;
            break;
        }

        s = Scheme::update(s);

    } while(true);

    detail::final_report(stream, msg, s.iteration, nfev_total, s.fx, s.dfx, dx);
    return std::make_pair(msg, x);
}

} // ns ook
#endif