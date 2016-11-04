#ifndef LBFGS_REPORT_HPP_
#define LBFGS_REPORT_HPP_

#include <array>
#include <cstdio>

#include "lbfgs_options.hpp"

namespace ook
{
namespace detail
{

/// \brief Prints a vector using the original LBFGS formatting.
/// \tparam A Vector type.
template <typename T>
void
print_vector(const T& x)
{
    int n = x.size();
    for (int i = 0; i < n; ++i)
    {
        printf("   %.3E", x[i]);
        if (i % 6 == 0)
            printf("\n");
    }
    printf("\n");
}

/// \brief Holder struct for some format strings.
namespace fmts
{
const constexpr char stars[] =
    "*************************************************\n";
const constexpr char header[] =
    "\n   I   NFN    FUNC        GNORM       STEPLENGTH\n\n";
}

/// \brief Print the output of required at the start of the iteration.
/// \tparam T Numeric type nested in lbfgs_options.
/// \tparam X Vector type.
/// \param opts Contains iprint variable.
/// \param gnorm The norm of the gradient.
/// \param x The current point.
/// \param f The value of the objective function at x.
/// \param g The gradient of the objective function at x.
template <typename T, typename X>
void
initial_report(
    const lbfgs_options<T>& opts, T gnorm, const X& x, T f, const X& g)
{
    if (opts.iprint[0] < 0)
        return;
    printf(fmts::stars);
    printf("  N=%5lu   NUMBER OF CORRECTIONS=%2d\n"
           "       INITIAL VALUES\n",
           x.size(),
           opts.m);
    printf(" F=  %.3E   GNORM=  %.3E\n", f, gnorm);
    if (opts.iprint[1] >= 1)
    {
        printf(" VECTOR X= \n");
        print_vector(x);
        printf("\n");
        printf(" GRADIENT VECTOR G= \n");
        print_vector(g);
    }
    printf(fmts::stars);
    printf(fmts::header);
}

/// \brief Print output according to the iprint values contained
/// in opts.options.
/// \tparam T Numeric type nested in lbfgs_options.
/// \tparam X Vector type.
/// \param opts lbfgs_options structure.
/// \param iter Current iteration.
/// \param nfun Number of function evaluations.
/// \param gnorm Norm of the gradient.
/// \param x Current point.
/// \param f Current function value.
/// \param g Current gradient.
/// \param stp Current Step.
/// \param finish Indicate if the algorithm has finished.
template <typename T, typename X>
void
report(const lbfgs_options<T>& opts,
       int iter,
       int nfun,
       T gnorm,
       const X& x,
       T f,
       const X& g,
       T stp,
       bool finish)
{
    auto iprint = opts.iprint;
    if (iprint[0] < 0)
        return;

    const char iterate[] = "%4d %4d    %10.3E  %10.3E  %10.3E\n";

    if (iprint[0] == 0 && (iter != 1 && !finish))
        return;

    if (iprint[0] != 0)
    {
        if ((iter - 1) % iprint[0] == 0 || finish)
        {
            if (iprint[1] > 1 && iter > 1)
            {
                printf(fmts::header);
            }
            printf(iterate, iter, nfun, f, gnorm, stp);
        }
    }
    else
    {
        if (iprint[1] > 1 && finish)
        {
            printf(fmts::header);
        }
        printf(iterate, iter, nfun, f, gnorm, stp);
    }

    if (iprint[1] == 2 || iprint[1] == 3)
    {
        if (finish)
        {
            printf(" FINAL POINT X= \n");
        }
        else
        {
            printf(" VECTOR X= \n");
        }
        print_vector(x);
        if (iprint[1] == 3)
        {
            printf(" GRADIENT VECTOR G= \n");
            print_vector(g);
        }
    }
    if (finish)
    {
        printf("\n THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS."
               "\n IFLAG = 0\n");
    }
}
} // ns detail
} // ns lbfgs

#endif
