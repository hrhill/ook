#include <iostream>
#include <string>

#include <boost/numeric/ublas/vector.hpp>

#include "ook.hpp"

typedef boost::numeric::ublas::vector<double> vector_type;

struct rosenbrock
{
    std::tuple<double, vector_type>
    operator()(const vector_type& x) const
    {
        const double x1 = x(0);
        const double x2 = x(1);
        const double f1 = 10 * (x2 - x1 * x1);
        const double f2 = 1- x1;
        const double f = f1 * f1 + f2 * f2;

        vector_type df(2, 1.0);
        df(0) = - 2 * f1 * 20 * x1 - 2 * f2;
        df(1) = 2 * f1 * 10;

        return std::make_tuple(f, df);
    }

    static const int n = 2;
};

template <typename State>
struct observer
{
    observer(const std::string& name)
    :
        name(name)
    {}

    void operator()(const State& s)
    {
        states.push_back(s);
    }

    friend
    std::ostream& operator<<(std::ostream& out, const observer& o)
    {
        out << o.name << o.states.back() << "\n";
        return out;
    }

    std::string name;
    std::vector<State> states;
};


int main(int argc, char** argv)
{
    using namespace ook;
    options<double> opts;

    vector_type x(2);
    x(0) = -1.2;
    x(1) = 1.0;
    typedef nonlinear_cg_impl<vector_type, beta::fr> fr_scheme;
    typedef line_search_method<fr_scheme>::state_type fr_state_type;
    observer<fr_state_type> fr_obs("fletcher-reeves");
    auto fr_soln = fletcher_reeves(rosenbrock(), x, opts, fr_obs);

    x(0) = -1.2;
    x(1) = 1.0;
    typedef nonlinear_cg_impl<vector_type, beta::pr> pr_scheme;
    typedef line_search_method<pr_scheme>::state_type pr_state_type;
    observer<pr_state_type> pr_obs("polak-ribiere");
    auto pr_soln = polak_ribiere(rosenbrock(), x, opts, pr_obs);

    x(0) = -1.2;
    x(1) = 1.0;
    typedef nonlinear_cg_impl<vector_type, beta::hs> hs_scheme;
    typedef line_search_method<hs_scheme>::state_type hs_state_type;
    observer<hs_state_type> hs_obs("hestenes-steifel");
    auto hs_soln = hestenes_steifel(rosenbrock(), x, opts, hs_obs);

    x(0) = -1.2;
    x(1) = 1.0;
    typedef nonlinear_cg_impl<vector_type, beta::dy> dy_scheme;
    typedef line_search_method<dy_scheme>::state_type dy_state_type;
    observer<dy_state_type> dy_obs("dai-yuan");
    auto dy_soln = dai_yuan(rosenbrock(), x, opts, dy_obs);

    std::cout << fr_obs << std::endl;
    std::cout << pr_obs << std::endl;
    std::cout << hs_obs << std::endl;
    std::cout << dy_obs << std::endl;
}