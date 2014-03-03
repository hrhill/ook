#ifndef OOK_LINE_SEARCH_MORE_THUENTE_SEARCH_H_
#define OOK_LINE_SEARCH_MORE_THUENTE_SEARCH_H_

#include <string>
#include <cmath>

#include "ook/line_search/more_thuente/step.h"
#include "ook/line_search/conditions.h"

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

//   stp is the current estimate of a satisfactory step.
//   f is the value of the function at stp.
//   g is the derivative of the function at stp.
namespace ook{
namespace line_search{
namespace detail{
struct state
{
    bool brackt;
    int stage;
    double gx;
    double gy;
    double fx;
    double fy;
    double stx;
    double sty;
    double stmin;
    double stmax;
    double width;
    double width1;
};

template <typename T, typename Options>
T
search(const T& istp, const T& f, const T& g, const T& f0, const T& g0, const Options& opts, state& s)
{
    T stp(istp);

    T gx = s.gx;
    T gy = s.gy;
    T fx = s.fx;
    T fy = s.fy;
    T stx = s.stx;
    T sty = s.sty;
    T stmin = s.stmin;
    T stmax = s.stmax;
    T width = s.width;
    T width1 = s.width1;

    const T xtrapl(1.1);
    const T xtrapu(4.0);

    const T gtest = opts.ftol * g0;
    const bool sufficient_decrease = sufficient_decrease_condition(f, f0, opts.ftol, stp, g0);
    // If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
    // algorithm enters the second stage.
    if (s.stage == 1 && sufficient_decrease && g >= 0.0) {
        s.stage = 2;
    }

    // A modified function is used to predict the step during the
    // first stage if a lower function value has been obtained but
    // the decrease is not sufficient.
    if (s.stage == 1 && f <= fx && !sufficient_decrease) {
        //    Define the modified function and derivative values.
        T fm = f - stp * gtest;
        T fxm = fx - stx * gtest;
        T fym = fy - sty * gtest;
        T gm = g - gtest;
        T gxm = gx - gtest;
        T gym = gy - gtest;
        // Call dcstep to update stx, sty, and to compute the new step.
        stp = step(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, s.brackt, stmin, stmax);
        // Reset the function and derivative values for f.
        fx = fxm + stx * gtest;
        fy = fym + sty * gtest;
        gx = gxm + gtest;
        gy = gym + gtest;
    } else {
        // Call dcstep to update stx, sty, and to compute the new step.
        stp = step(stx, fx, gx, sty, fy, gy, stp, f, g, s.brackt, stmin, stmax);
    }
    // Decide if a bisection step is needed.
    if (s.brackt) {
        if (fabs(sty - stx) >= 0.66 * width1) {
            stp = stx + 0.5 * (sty - stx);
        }
        width1 = width;
        width = fabs(sty - stx);
    }
    // Set the minimum and maximum steps allowed for stp.
    if (s.brackt) {
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
    if ((s.brackt && (stp <= stmin || stp >= stmax)) ||
        (s.brackt && (stmax - stmin <= opts.xtol * stmax)))
    {
        stp = stx;
    }

    // Save state
    s.gx = gx;
    s.gy = gy;
    s.fx = fx;
    s.fy = fy;
    s.stx = stx;
    s.sty = sty;
    s.stmin = stmin;
    s.stmax = stmax;
    s.width = width;
    s.width1 = width1;

    return stp;
}
/*
template <typename T>
struct search{
    search(const T& gx, const T& gy,
           const T& fx, const T& fy,
           const T& stx, const T& sty,
           const T& stmin, const T& stmax,
           const T& width, const T& width1)
    :
        brackt_(false),
        stage_(1),
        gx_(gx),
        gy_(gy),
        fx_(fx),
        fy_(fy),
        stx_(stx),
        sty_(sty),
        stmin_(stmin),
        stmax_(stmax),
        width_(width),
        width1_(width1)
    {}



private:
    bool brackt_;
    int stage_;
    T gx_;
    T gy_;
    T fx_;
    T fy_;
    T stx_;
    T sty_;
    T stmin_;
    T stmax_;
    T width_;
    T width1_;
};
*/
} // detail
} // ns line_search
} // ns ook

#endif
