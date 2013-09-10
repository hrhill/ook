/// \brief Implementation of the functions described in 
/// "Testing unconstrained optimization software", 
/// More JJ, Garbow BS, Hillstrom KE
#ifndef MORE_GARBOW_HILLSTROM_TEST_FUNCTIONS_H_
#define MORE_GARBOW_HILLSTROM_TEST_FUNCTIONS_H_

#include <tuple>
#include <vector>

template <typename Vector>
struct rosenbrock{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        const real_type f1 = 10.0 * (x(1) - x(0) * x(0));
        const real_type f2 = 1.0 - x(0);
        const real_type f = f1 * f1 + f2 * f2;
        vector_type df(2);
        df(0) = -4.0 * 10.0 * f1 * x(0) - 2.0 * f2;
        df(1) = -2.0 * 10.0 * f1;
        return std::make_pair(f, df); 
    }

    static const int n = 2;
    static const int m = 2;
    static real_type f_min;
    static std::vector<real_type> minima;
    static std::vector<real_type> x0;    
};


template <typename Vector>
std::vector<typename Vector::value_type>
rosenbrock<Vector>::minima = {1.0, 1.0};

template <typename Vector>
std::vector<typename Vector::value_type>
rosenbrock<Vector>::x0 = {-1.2, 1.0};

template <typename Vector>
typename rosenbrock<Vector>::real_type
rosenbrock<Vector>::f_min = 0.0;

/*
template <typename Vector>
struct freudenstein_roth{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct powell_badly_scaled{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct brown_badly_scaled{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct beale{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct jenrich_sampson{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct helical_valley{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct bard{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct gaussian{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct meyer{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct gulf_rnd{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct box_3d{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct powell_singular{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct wood{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct kowalik_osborne{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct brown_dennis{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct osborne_1{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct biggs_exp6{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct osborne_2{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct watson{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct extended_rosenbrock{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct extended_powell_singular{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct penalty_i{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct penalty_ii{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct variable_dim{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct trigonometric{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct brown_almost_linear{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct discrete_boundary_value{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct discrete_integral_equation{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct broyden_tridiagonal{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct broyden_banded{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct linear_full_rank{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct linear_rank_1{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct linear_rank_1_with_0_cols_rows{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};

template <typename Vector>
struct chebyquad{

    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        return std::make_pair(0.0, x); 
    }
};
*/
#endif