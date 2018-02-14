#include <iostream>
#include <string>

#include <boost/numeric/ublas/vector.hpp>
#include <utility>

#include "ook.hpp"

struct rosenbrock
{
    std::tuple<double, ook::vector>
    operator()(const ook::vector& x) const
    {
        const double x1 = x[0];
        const double x2 = x[1];
        const double f1 = 10 * (x2 - x1 * x1);
        const double f2 = 1 - x1;
        const double f = f1 * f1 + f2 * f2;

        ook::vector df(2, 1.0);
        df[0] = -2 * f1 * 20 * x1 - 2 * f2;
        df[1] = 2 * f1 * 10;

        return std::make_tuple(f, df);
    }

    static const int n = 2;
};

template <typename State>
struct observer
{
    explicit observer(std::string  name) : name(std::move(name)) {}

    void
    operator()(const State& s)
    {
        states.push_back(s);
    }

    friend std::ostream&
    operator<<(std::ostream& out, const observer& o)
    {
        out << o.name << o.states.back() << "\n";
        return out;
    }

    std::string name;
    std::vector<State> states;
};

int
main()
{
    using namespace ook;
    options<double> opts;

    ook::vector x(2);
    x[0] = -1.2;
    x[1] = 1.0;
    using fr_scheme = nonlinear_cg_impl<beta::fr>;
    typedef line_search_method<fr_scheme, line_search::mcsrch>::state_type
        fr_state_type;
    observer<fr_state_type> fr_obs("fletcher-reeves");
    auto fr_soln = fletcher_reeves(rosenbrock(), x, opts, fr_obs);

    x[0] = -1.2;
    x[1] = 1.0;
    using pr_scheme = nonlinear_cg_impl<beta::pr>;
    typedef line_search_method<pr_scheme, line_search::mcsrch>::state_type
        pr_state_type;
    observer<pr_state_type> pr_obs("polak-ribiere");
    auto pr_soln = polak_ribiere(rosenbrock(), x, opts, pr_obs);

    x[0] = -1.2;
    x[1] = 1.0;
    using hs_scheme = nonlinear_cg_impl<beta::hs>;
    typedef line_search_method<hs_scheme, line_search::mcsrch>::state_type
        hs_state_type;
    observer<hs_state_type> hs_obs("hestenes-steifel");
    auto hs_soln = hestenes_steifel(rosenbrock(), x, opts, hs_obs);

    x[0] = -1.2;
    x[1] = 1.0;
    using dy_scheme = nonlinear_cg_impl<beta::dy>;
    typedef line_search_method<dy_scheme, line_search::mcsrch>::state_type
        dy_state_type;
    observer<dy_state_type> dy_obs("dai-yuan");
    auto dy_soln = dai_yuan(rosenbrock(), x, opts, dy_obs);

    std::cout << fr_obs << std::endl;
    std::cout << pr_obs << std::endl;
    std::cout << hs_obs << std::endl;
    std::cout << dy_obs << std::endl;
}
