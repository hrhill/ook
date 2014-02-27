#ifndef OOK_LINE_SEARCH_METHODS_LBFGS_H_
#define OOK_LINE_SEARCH_METHODS_LBFGS_H_

#include "ook/lbfgs/lbfgs.h"

namespace ook{

template <typename F, typename X, typename Options, typename Observer>
std::tuple<ook::message, X>
lbfgs(F objective_function, X x, const Options& opts, Observer& observer)
{
    const int n = x.size();
    int m = 5;

    X diag(n);
    X g(n);
    X w(n * (2 * m + 1) + 2 * m);

    int iprint[2] = {1, 0};

    bool diagco = false;
    double eps = 1e-5;
    int icall = 0;
    int iflag = 0;

    std::fill(diag.begin(), diag.end(), 1.0);

    while(true){
        double f = 0.;
        std::tie(f, g) = objective_function(x);

        detail::lbfgs(n, m, &x[0], f, &g[0], diagco, &diag[0], iprint, eps, opts.xtol, &w[0], iflag);
        if (iflag <= 0) {
            break;
        }
        ++icall;
        if (icall > opts.max_function_evaluations) {
            break;
        }
    }
    return std::make_tuple(ook::message::convergence, x);
}

} // ns ook

#endif
