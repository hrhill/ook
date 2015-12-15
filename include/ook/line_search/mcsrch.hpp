// Copyright 2013 Harry Hill
//
// This file is part of ook.
//
// ook is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// ook is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with ook.  If not, see <http://www.gnu.org/licenses/>.

#ifndef OOK_LINE_SEARCH_MORE_THUENTE_MCSRCH_HPP_
#define OOK_LINE_SEARCH_MORE_THUENTE_MCSRCH_HPP_

#include <limits>
#include "ook/message.hpp"

#include "./mcstep.hpp"

namespace ook{
namespace line_search{

/// \brief The purpose of mcsrch is to find a step which satisfies a sufficient
/// decrease condition and a curvature condition.
/// \detail At each stage the function updates an interval of uncertainty with
/// endpoints stx and sty. The interval of uncertainty is initially chosen so
/// that it contains a minimizer of the modified function
/// \f[
///    f(x + stp * s) - f(x) - ftol * stp * grad(x)'s
/// \f]
/// If a step is obtained for which the modified function has a non-positive
/// function value and non-negative derivative , then the interval of
/// uncertainty is chosen so that is contains a minimizer of f(x + stp * s)
///
/// The algorithm is designed to find a step which satisfies the sufficient
/// decrease condition,
/// \f[
/// f(x + stp s) <= f(x) + ftol stp  \nabla f(x)' s
/// \f]
/// and the curvature condition.
/// \f[
/// |  |\nabla f(x + stp * s)'s| <= gtol |\nabla f(x)'s)|.
/// \f]
/// If ftol is less gtol, and if, for example the function is bounded below,
/// then there is always a step which satisfies both conditions. If no step can
/// be found which satisfies both conditions, then the algorithm usually stops
/// when rounding errors prevent further progress. In this case, stp usually
/// only satisfies the sifficient decrease condition.
///
/// \tparam F Callable type, with signature std::tuple<T, T> (T);
/// \tparam T Numeric type supporting +,-,*, /, <, >, fabs
/// \param phi One dimensional function to minimize.
/// \param finit The function value at the initial point.
/// \param dginit The derivative at the initial point.
/// \param stp The initial step value.
/// \param opts lbfgs_options structure.
/// \return info, stp The values of info and the step.
///   info = 0  Improper input parameters.
///   info = 1  The sufficient decrease and the directional derivative
///             condition hold.
///   info = 2  Relative width of the interval of uncertainty is at most xtol.
///   info = 3  Number of calls to phi has reached maxfev.
///   info = 4  The step is at the lower bound stpmin.
///   info = 5  The step is at the lower bound stpmax.
///   info = 6  Rounding errors prevent further progress. There may not be a
///             step which satisfies the sufficient decrease and curvature
///             conditions. Tolerances may be too small.
struct mcsrch
{
    template <typename F, typename T, typename Options>
    std::tuple<ook::message, T, T, T>
    operator()(F phi, T finit, T dginit, T stp, const Options& opts) const
    {
        const T ftol = opts.ftol;
        const T xtol = std::numeric_limits<T>::epsilon();
        const T p5 = .5;
        const T p66 = .66;
        const T xtrapf = 4.;

        const T stpmin = 1e-20;
        const T stpmax = 1e20;

        int infoc = 1;
        int nfev = 0;
        bool brackt = false;
        bool stage1 = true;
        const T dgtest = ftol * dginit;
        T width = stpmax - stpmin;
        T width1 = width / p5;

        // The values stx, fx, dgx are step, function and derivative at the best
        // step.
        // The values sty, fy, dgy are the step, function and derivative at the
        // other endpoint of the interval of uncertainty.
        // The values stp, f, dg are the step, function and derivative at the
        // current step.
        T stx = 0.0;
        T fx = finit;
        T dgx = dginit;
        T sty = 0.0;
        T fy = finit;
        T dgy = dginit;

        while(true)
        {
            T stmin = stx;
            T stmax = stp + xtrapf * (stp - stx);

            if(brackt){
                // Set max and min step to the present interval of uncertainty.
                stmin = std::min(stx, sty);
                stmax = std::max(stx, sty);
            }

            // Clamp the step to be between maximum and minimum values.
            stp = std::max(stp, stpmin);
            stp = std::min(stp, stpmax);

            // If an unusual termination occues, set stp to be the lowest point
            // obtained so far.
            if((brackt && (stp <= stmin || stp >= stmax)) || infoc == 0 ||
                (brackt && stmax - stmin <= xtol * stmax))
            {
                stp = stx;
            }
            // Evaluate the function and derivative.
            T f, dg;
            std::tie(f, dg) = phi(stp);
            const T ftest1 = finit + stp * dgtest;

            // Test for convergence.
            ook::message msg = ook::message::null;
            if((brackt && (stp <= stmin || stp >= stmax)) || infoc == 0){
                msg = ook::message::warning_rounding_error_prevents_progress;
            }
            if(stp == stpmax && f <= ftest1 && dg <= dgtest){
                msg = ook::message::warning_stp_eq_stpmax;
            }
            if(stp == stpmin && (f > ftest1 || dg >= dgtest)){
                msg = ook::message::warning_stp_eq_stpmin;
            }
            if(nfev >= opts.maxfev){
                msg = ook::message::warning_max_line_search_attempts_reached;
            }
            if(brackt && stmax - stmin <= xtol * stmax){
                msg = ook::message::warning_xtol_satisfied;
            }
            if(f <= ftest1 && fabs(dg) <= opts.gtol * (-dginit)){
                msg = ook::message::convergence;
            }
            // Check for termination.
            if(msg != ook::message::null){
                return std::make_tuple(msg, stp, f, dg);
            }

            // In the first state, seek a step for which the modified function has
            // a non-positive value and non-negative derivative.
            if(stage1 && f <= ftest1 && dg >= std::min(ftol, opts.gtol) * dginit){
                stage1 = false;
            }
            // A modified function is used to predict the step only if we have not
            // obtained a step for which the modified function has a non-positive
            // function value and non-negative derivative, and if a lower function
            // value has been obtained but the decrease is not sufficient.
            if(stage1 && f <= fx && f > ftest1){
                /// Define the modified function and derivative values.
                T fm = f - stp * dgtest;
                T fxm = fx - stx * dgtest;
                T fym = fy - sty * dgtest;
                T dgm = dg - dgtest;
                T dgxm = dgx - dgtest;
                T dgym = dgy - dgtest;

                // Compute new step and update the interval of uncertainty.
                infoc = mcstep(stx, fxm, dgxm,
                                sty, fym, dgym,
                                stp, fm, dgm, brackt, stmin, stmax);

                // Reset function and gradient values.
                fx = fxm + stx * dgtest;
                fy = fym + sty * dgtest;
                dgx = dgxm + dgtest;
                dgy = dgym + dgtest;
            } else {
                // Compute new step and update the interval of uncertainty.
                infoc = mcstep(stx, fx, dgx,
                                sty, fy, dgy,
                                stp, f, dg, brackt, stmin, stmax);
            }
            // Force a sufficient decrease in the size od the interval of
            // uncertainty.
            if(brackt) {
                if(fabs(sty - stx) >= p66 * width1)
                {
                    stp = stx + p5 * (sty - stx);
                }
                width1 = width;
                width = fabs(sty - stx);
            }
        }
    }
};

} // ns line_search
} // ns ook

#endif
