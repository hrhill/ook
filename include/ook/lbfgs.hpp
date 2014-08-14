#ifndef LBFGS_HPP_
#define LBFGS_HPP_

#include <cmath>
#include <array>
#include <cassert>
#include <algorithm>

extern "C" {
#include <cblas.h>
}
#include "ook/line_search/mmt/mcsrch.hpp"
#include "ook/lbfgs/report.hpp"

namespace ook{

/// \brief Limited memory BFGS method for large scale optimisation,
/// Jorge Nocedal, July 1990.
/// \detail This subroutine solves the unconstrained minimization problem
/// \f[
///                  \mbox{min} f(x),\quad    x \in \mathbb{R}^n,
/// \f]
/// using the limited memory BFGS method. The routine is especially effective
/// on problems involving a large number of variables. In a typical iteration
/// of this method an approximation Hk to the inverse of the Hessian is
/// obtained by applying M BFGS updates to a diagonal matrix Hk0, using
/// information from the previous M steps.
///
/// The user specifies the number M, which determines the amount of storage
/// required by the routine. The user may also provide the diagonal matrices
/// Hk0 if not satisfied with the default choice. The algorithm is described in
///     "On the limited memory BFGS method for large scale optimization",
///     by D. Liu and J. Nocedal,
/// Mathematical Programming B 45 (1989) 503-528.
///
/// The steplength is determined at each iteration by means of the line search
/// routine mcsrch, which is a slight modification of the routine csrch written
/// by More' and Thuente.
///
/// \tparam F Type of objective function.
/// \tparam D Type of function that computes the diagonal.
/// \tparam X Vector type.
/// \tparam T Numeric type, nested in the options structure.
/// \param obj_f The objective function to be minimized.
/// \param diag_f A function returning the diagonal matrix Hk0.
/// \param x The initial point.
/// \param opts A structure containing the required options.
/// \return A flag together with the current best solution.
template <typename F, typename D, typename X, typename T>
std::tuple<int, X>
lbfgs(F obj_f, D diag_f, X x, const lbfgs_options<T>& opts)
{
    const int n = x.size();

    // Initialize function and gradient.
    double f = 0;
    X g(n);
    std::tie(f, g) = obj_f(x);

    X diag(n);
    if(opts.diagco){
        diag = diag_f(x);
    }else{
        std::fill(diag.begin(), diag.end(), 1.0);
    }

    // The work vector is divided as follws:
    // - The first n locations are used to store the gradient and other
    // temporary information
    // - locations (N+1)...(N+M) store the scalars rho.
    // - locations (N+M+1)...(N+2M) store the numbers alpha used in the formula
    //   that computes H * G.
    // - locations (N+2M+1)...(N+2M+NM) store the last m steps.
    // - locations (N+2M+NM+1)...(N+2M+2NM) store the last m gradient
    //   differences.
    // The search steps and gradient differences are stored in a circular
    // order controlled bt the parameter point.
    const int m = opts.m;
    const int wa_size = n * (2 * m + 1) + 2 * m;
    X w(wa_size);

    int ispt = n + (m << 1);
    int iypt = ispt + n * m;
    for(int i = 0; i < n; ++i) {
        w[ispt + i] = -g[i] * diag[i];
    }
    double gnorm = sqrt(cblas_ddot(n, &g[0], 1, &g[0], 1));
    const double stp1 = 1.0 / gnorm;

    detail::initial_report(opts, gnorm, x, f, g);

    int iter = 0;
    int nfun = 1;
    int point = 0;
    int npt = 0;

    while(true){
        ++iter;
        int info = 0;
        int bound = iter - 1;

        if(iter != 1){
            if(iter > m){
                bound = m;
            }

            const double ys = cblas_ddot(n, &w[iypt + npt], 1,
                                            &w[ispt + npt], 1);
            if (opts.diagco){
                diag = diag_f(x);
            } else {
                const double yy = cblas_ddot(n, &w[iypt + npt], 1,
                                                &w[iypt + npt], 1);
                std::fill(diag.begin(), diag.end(), ys/yy);
            }
            // Compute -H*G using the formula given in : Nocedal, J. 1980,
            // "Updating quasi-Newton matrices with limited storage",
            // Mathematics of Computation, Vol.24, No.151, pp. 773-782.
            int cp = point;
            if(point == 0){
                cp = m;
            }
            w[n + cp - 1] = 1.0 / ys;
            for(int i = 0; i < n; ++i){
                w[i] = -g[i];
            }
            cp = point;
            int inmc = 0;
            for(int i = 0; i < bound; ++i){
                --cp;
                if(cp == -1) {
                    cp = m - 1;
                }
                const double sq = cblas_ddot(n, &w[ispt + cp * n], 1,
                                                &w[0], 1);
                inmc = n + m + cp + 1;
                int iycn = iypt + cp * n;
                w[inmc - 1] = w[n + cp] * sq;
                double mwinmc = -w[inmc - 1];
                cblas_daxpy(n, mwinmc, &w[iycn], 1, &w[0], 1);
            }

            for(int i = 0; i < n; ++i){
                w[i] = diag[i] * w[i];
            }

            for(int i = 0; i < bound; ++i){
                const double yr = cblas_ddot(n, &w[iypt + cp * n], 1,
                                                &w[0], 1);
                double beta = w[n + cp] * yr;
                inmc = n + m + cp + 1;
                beta = w[inmc - 1] - beta;
                const int iscn = ispt + cp * n;
                cblas_daxpy(n, beta, &w[iscn], 1, &w[0], 1);
                ++cp;
                cp %= m;
            }
            // Store the new search direction.
            std::copy(w.begin(), w.begin() + n, w.begin() + ispt + point * n);
        }
        // Obtain the 1-dimensional minimizer of the line search function.
        double stp = (iter == 1) ? stp1 : 1.0;

        std::copy(g.begin(), g.end(), w.begin());
        X s(w.begin() + ispt + point * n, w.begin() + ispt + point * n + n);

        double dginit = cblas_ddot(n, &g[0], 1, &s[0], 1);
        if(dginit >= 0.0)
            throw std::runtime_error(
                        "The search direction is not a descent direction.");

        X x0(x);
        auto phi = [&nfun, &f, &x0, &x, &g, obj_f, &s](double stp)
        {
            ++nfun;
            int n = x.size();
            for (int j = 0; j < n; ++j){
                x[j] = x0[j] + stp * s[j];
            }
            std::tie(f, g) = obj_f(x);
            double dg = cblas_ddot(n, &g[0], 1, &s[0], 1);
            return std::make_tuple(f, dg);
        };
        std::tie(info, stp) = line_search::mmt::mcsrch(phi, f, dginit, stp, opts);

        if(info != 1){
            printf(" IFLAG= -1\n LINE SEARCH FAILED."
             " SEE DOCUMENTATION OF ROUTINE MCSRCH\n"
             " ERROR RETURN OF LINE SEARCH: INFO= %2d\n"
             " POSSIBLE CAUSES: FUNCTION OR GRADIENT ARE INCORRECT\n"
             " OR INCORRECT TOLERANCES\n", info);
            return std::make_tuple(-1, x);
        }

        // Compute the new step and gradient change.
        npt = point * n;
        for(int i = 0; i < n; ++i) {
            w[ispt + npt + i] *= stp;
            w[iypt + npt + i] = g[i] - w[i];
        }
        ++point;
        if(point == m) {
            point = 0;
        }

        // Termination test
        gnorm = sqrt(cblas_ddot(n, &g[0], 1, &g[0], 1));
        const double xnorm = std::max(1.0,
                                    sqrt(cblas_ddot(n, &x[0], 1, &x[0], 1)));
        bool finish = (gnorm / xnorm <= opts.eps) || (iter >= opts.maxicall);
        detail::report(opts, iter, nfun, gnorm, x, f, g, stp, finish);

        if(finish){
            break;
        }
    }
    return std::make_tuple(0, x);
}

template <typename F, typename X, typename T>
std::tuple<int, X>
lbfgs(F obj_f, X x, const lbfgs_options<T>& opts)
{
    auto diag_f = [](const X& x){
        return X();
    };
    return lbfgs(obj_f, diag_f, x, opts);
}


} // ns ook

#endif
