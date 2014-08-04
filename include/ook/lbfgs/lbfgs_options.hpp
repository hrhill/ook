#ifndef LBFGS_LBFGS_OPTIONS_HPP_
#define LBFGS_LBFGS_OPTIONS_HPP_

#include <array>
#include <cassert>

namespace ook
{
/// \brief User controlled options used by lbfgs and the line search algorithm.
template <typename T>
struct lbfgs_options
{
    /// \brief Explicit constructor.
    lbfgs_options(int m,
                  T eps,
                  bool diagco,
                  std::array<int, 2> iprint,
                  T ftol,
                  T gtol,
                  int maxicall,
                  int maxfev)
    :
        m(m),
        eps(eps),
        diagco(diagco),
        iprint(iprint),
        ftol(ftol),
        gtol(gtol),
        maxicall(maxicall),
        maxfev(maxfev)
    {
        assert(m > 0);
        assert(maxfev > 0);
        assert(gtol > 1e-04);
    }

    /// \brief The number of corrections used in the BFGS update. Values less
    /// than 3 are not recommended; large values of m will result in excessive
    /// computing time. 3 <= m <= 7 is recommended. Typical choice is m = 5.
    int m;

    /// \brief Determines the accuracy with which the solution is to be found.
    /// The algorithm terminates when \f$ \|G\| < \epsilon \mbox{max}(1, \|X\|)\f$.
    /// A typical choice is 1e-05.
    T eps;

    /// \brief Set to true if the user wishes to provide the diagonal matrix
    /// Hk0 at each iteration. Otherwise it should be set to false, in which
    /// case LBFGS will use a default value described in the algorithm. If
    /// diagco is set to true the routine will call the user supplied function.
    bool diagco;

    /// \brief An integer array of length two which must be set by the user.
    ///
    /// iprint[0] specifies the frequency of the output:
    ///     iprint[0] < 0 : no output is generated,
    ///     iprint[0] = 0 : output only at first and last iteration,
    ///     iprint[0] > 0 : output every iprint[0] iterations.
    ///
    /// iprint[1] specifies the type of output generated:
    ///     iprint[1] = 0 : iteration count, number of function evaluations,
    ///                     function value, norm of the gradient, and
    ///                     steplength,
    ///     iprint[1] = 1 : same as iprint[1]=0, plus vector of variables and
    ///                     gradient vector at the initial point,
    ///     iprint[1] = 2 : same as iprint[1]=1, plus vector of variables,
    ///     iprint[1] = 3 : same as iprint[1]=2, plus gradient vector.
    /// Typical value is {1, 0}.
    std::array<int, 2> iprint;

    /// \brief The constant in the sufficient decrease condition. Typical value
    /// 1e-04.
    T ftol;

    /// \brief The constant in the curvature condition, typical value 0.9. Must
    /// be greater than 1e-04.
    T gtol;

    /// \brief The maximum number of iterations. Typical value 2000.
    int maxicall;

    /// \brief Maximum number of function evaluations used in the line search.
    /// Typical value 20.
    int maxfev;
};

} // ns lbfgs

#endif
