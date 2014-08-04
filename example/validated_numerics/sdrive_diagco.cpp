#include <vector>
#include <tuple>

#include "ook/lbfgs.hpp"

std::tuple<double, std::vector<double>>
obj_fun(const std::vector<double>& x)
{
    const int n = x.size();
    std::vector<double> g(n);
    double f = 0.;
    for (int j = 0; j < n; j += 2){
        double t1 = 1.0 - x[j];
        double t2 = (x[j + 1] - x[j] * x[j]) * 10.;
        g[j + 1] = t2 * 20.0;
        g[j] = (x[j] * g[j + 1] + t1) * -2.0;
        f += t1 * t1 + t2 * t2;
    }
    return std::make_tuple(f, g);
}

std::vector<double>
diag_fun(const std::vector<double>& x)
{
    const int n = x.size();
    std::vector<double> diag(n);
    for (int j = 0; j < n; ++j)
    {
        diag[j] = 1.0;
    }
    return diag;
}

int main()
{
    ook::lbfgs_options<double>
    opts(5, 1e-05, true, {1, 0}, 1e-04, 0.9, 2000, 20);

    // Rosenbrock problem, only valid for even n.
    int n = 100;
    std::vector<double> x(n);
    for (int j = 0; j < n; j += 2){
        x[j] = -1.2;
        x[j + 1] = 1.0;
    }

    auto soln = ook::lbfgs(obj_fun, diag_fun, x, opts);
}
