#ifndef OOK_LINE_SEARCH_METHOD_H_
#define OOK_LINE_SEARCH_METHOD_H_

#include <boost/numeric/ublas/matrix.hpp>

namespace ook{

namespace detail{

template <typename X>
struct 
state{
    typedef X vector_type;
    typedef typename X::value_type value_type;
    
    typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> matrix_type;
    state(const int n)
    :
         dfx(n), dfx0(n), p(n), H(boost::numeric::ublas::identity_matrix<double>(n)), a(1), beta(0), iteration(0)
    {}


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

template <typename X>
void
report(int iteration, int nfev_total, int nfev, double a, double fx, const X& dfx, const X& dx)
{
    std::cout << std::setw(8) << iteration
                 << std::setw(16) << nfev_total
                 << std::setw(8) << nfev 
                 << std::scientific 
                 << std::setw(16) << a
                 << std::setw(16) << fx
                 << std::setw(16) << ook::norm_infinity(dfx)
                 << std::setw(16) << ook::norm_infinity(dx) << std::endl;  
}

template <typename X>
void
final_report(int nfev_total,double fx, const X& dfx, const X& dx)
{
    std::cout << std::setw(8) << std::scientific 
                 << std::setw(16) << nfev_total
                 << std::setw(16) << fx
                 << std::setw(16) << ook::norm_infinity(dfx)
                 << std::setw(16) << ook::norm_infinity(dx) << std::endl;  
}


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

template <typename Scheme, typename F, typename X, typename Options>
std::tuple<ook::state_value, X>
line_search_method(F objective_function, X x, const Options& opts)
{
    typedef typename X::value_type real_type;
    typedef typename Scheme::state_type state_type;

    typedef decltype(objective_function(x)) result_type;
    typedef detail::function_caller<F, X, state_type, std::tuple_size<result_type>::value> fcaller_type;

    state_type s = Scheme::initialise(objective_function, x);

    s.iteration = 0;
    uint nfev_total = 0;
    ook::state_value value;

    do {
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
        std::tie(value, s.a) = ook::line_search::more_thuente(phi, s.fx, dfx_dot_p, s.a, opts);

        X dx(s.a * p);
        x += dx;
        nfev_total += nfev;
        ++s.iteration;

        //detail::report(s.iteration, nfev_total, nfev, s.a, s.fx, s.dfx, dx);
        if (ook::norm_infinity(s.dfx) < 1e-08){
            value = ook::state_value::convergence;
            detail::final_report(nfev_total, s.fx, s.dfx, dx);
            break;
        }

        s = Scheme::update(s);
    } while(true);

    return std::make_pair(value, x);
}

} // ns ook
#endif