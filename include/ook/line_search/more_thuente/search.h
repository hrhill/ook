#ifndef OOK_LINE_SEARCH_MORE_THUENTE_SEARCH_H_
#define OOK_LINE_SEARCH_MORE_THUENTE_SEARCH_H_

#include <tuple>
#include <cmath>

#include "ook/line_search/more_thuente/step.h"
#include "ook/line_search/conditions.h"
#include "ook/message.h"

// This subroutine finds a step that satisfies a sufficient
// decrease condition and a curvature condition.

// Each call of the subroutine updates an interval with
// endpoints stx and sty. The interval is initially chosen
// so that it contains a minimizer of the modified function

//       psi(stp) = f(stp) - f(0) - ftolstp*f'(0).

// If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
// interval is chosen so that it contains a minimizer of f.

// The algorithm is designed to find a step that satisfies
// the sufficient decrease condition

//       f(stp) <= f(0) + ftolstp*f'(0),

// and the curvature condition

//       fabs(f'(stp)) <= gtol*fabs(f'(0)).

// If ftol is less than gtol and if, for example, the function
// is bounded below, then there is always a step which satisfies
// both conditions.

// If no step can be found that satisfies both conditions, then
// the algorithm stops with a warning. In this case stp only
// satisfies the sufficient decrease condition.
namespace ook{
namespace line_search{
namespace detail{

template <typename T>
struct search{

    search(T st0_, T f0_, T g0_, T sty_, T fy_, T gy_, T stmin_, T stmax_, T width_, T width1_)
    :
        bracket(false),
        stage(1),
        f0(f0_),
        g0(g0_),
        stx(st0_),
        fx(f0_),
        gx(g0_),
        sty(sty_),
        fy(fy_),
        gy(gy_),
        stmin(stmin_),
        stmax(stmax_),
        width(width_),
        width1(width1_)
    {}

    template <typename Options>
    std::tuple<message, T>
    operator()(T stp, T f, T g, const Options& opts)
    {
        const T xtrapl(1.1);
        const T xtrapu(4.0);

        const T gtest = opts.ftol * g0;
        // Test for convergence.
        const bool curvature = curvature_condition(g, opts.gtol, g0);
        const bool sufficient_decrease = sufficient_decrease_condition(f, f0, opts.ftol, stp, g0);

        if (sufficient_decrease && curvature)
            return std::make_tuple(message::convergence, stp);
        // Test for warnings.
        if (bracket && (stp <= stmin || stp >= stmax))
            return std::make_tuple(message::warning_rounding_error_prevents_progress, stp);
        if (bracket && stmax - stmin <= opts.xtol * stmax)
            return std::make_tuple(message::warning_xtol_satisfied, stp);
        if (stp == opts.stpmax && sufficient_decrease && curvature)
            return std::make_tuple(message::warning_stp_eq_stpmax, stp);
        if (stp == opts.stpmin && (!sufficient_decrease || !curvature))
            return std::make_tuple(message::warning_stp_eq_stpmin, stp);

        // If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
        // algorithm enters the second stage.
        if (stage == 1 && sufficient_decrease && g >= 0.0) {
            stage = 2;
        }

        // A modified function is used to predict the step during the
        // first stage if a lower function value has been obtained but
        // the decrease is not sufficient.
        if (stage == 1 && f <= fx && !sufficient_decrease) {
            //    Define the modified function and derivative values.
            T fm = f - stp * gtest;
            T fxm = fx - stx * gtest;
            T fym = fy - sty * gtest;
            T gm = g - gtest;
            T gxm = gx - gtest;
            T gym = gy - gtest;
            // Call dcstep to update stx, sty, and to compute the new step.
            stp = step(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, bracket, stmin, stmax);
            // Reset the function and derivative values for f.
            fx = fxm + stx * gtest;
            fy = fym + sty * gtest;
            gx = gxm + gtest;
            gy = gym + gtest;
        } else {
            // Call dcstep to update stx, sty, and to compute the new step.
            stp = step(stx, fx, gx, sty, fy, gy, stp, f, g, bracket, stmin, stmax);
        }
        // Decide if a bisection step is needed.
        if (bracket) {
            if (fabs(sty - stx) >= 0.66 * width1) {
                stp = stx + 0.5 * (sty - stx);
            }
            width1 = width;
            width = fabs(sty - stx);
        }
        // Set the minimum and maximum steps allowed for stp.
        if (bracket) {
            stmin = std::min(stx, sty);
            stmax = std::max(stx, sty);
        } else {
            stmin = stp + xtrapl * (stp - stx);
            stmax = stp + xtrapu * (stp - stx);
        }
        // Force the step to be within the bounds stpmax and stpmin.
        stp = std::max(stp, opts.stpmin);
        stp = std::min(stp, opts.stpmax);
        // If further progress is not possible, let stp be the best
        // point obtained during the search.
        if ((bracket && (stp <= stmin || stp >= stmax)) ||
            (bracket && (stmax - stmin <= opts.xtol * stmax)))
        {
            stp = stx;
        }
        return std::make_tuple(message::update, stp);
    }

private:
    bool bracket;
    int stage;

    T f0;
    T g0;

    T stx;
    T fx;
    T gx;

    T sty;
    T fy;
    T gy;

    T stmin;
    T stmax;
    T width;
    T width1;
};

} // detail
} // ns line_search
} // ns ook

#endif
