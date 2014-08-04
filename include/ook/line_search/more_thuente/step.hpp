// Copyright 2013 Harry Hill
//
// This file is part of ook.
//
// ook is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ook is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with ook.  If not, see <http://www.gnu.org/licenses/>.

#ifndef OOK_LINE_SEARCH_MORE_THUENTE_STEP_HPP_
#define OOK_LINE_SEARCH_MORE_THUENTE_STEP_HPP_

#include <cmath>
#include <algorithm>

namespace ook{
namespace line_search{
namespace detail{

/// \brief Sign function
template <typename T>
T
signum(T x)
{
    return copysign(1.0, x);
}

/// \brief Calculate the secant step.
template <typename T>
T
secant_step(T x, T dx, T y, T dy)
{
    return y + dy / (dy - dx) * (x - y);
}

/// \brief Calculate the value which minimizes the quadratic
/// interpolant to the inputs.
template <typename T>
T
quadratic_step(T x, T fx, T dx, T y, T fy)
{
    return  x + dx / ((fx - fy) / (y - x) + dx) / 2.0 * (y - x);
}

/// \brief Calculate the value which minimizes the cubic
/// interpolant to the inputs.
template <typename T>
T
cubic_step(T x, T fx, T dx, T y, T fy, T dy)
{
    const T d1 = dx + dy - 3 * (fx - fy)/(x - y);
    const T s = std::max({fabs(d1), fabs(dx), fabs(dy)});
    const T d2 = signum(y - x) * s * sqrt(pow(d1 / s, 2) - (dx / s) * (dy / s));
    return y - (y - x) * (dy + d2 - d1)/(dy - dx + 2 * d2);
}

/// \brief Return if a is closer to c than b.
template <typename T>
bool
is_closer(T a, T b, T c)
{
    return fabs(a - c) < fabs(b - c);
}

// This function computes a safeguarded step for a search
// procedure and updates an interval that contains a step that
// satisfies a sufficient decrease and a curvature condition.
//
// The parameter stx contains the step with the least function
// value. If bracket == true then a minimizer has
// been bracketed in an interval with endpoints stx and sty.
// The parameter stp contains the current step.
//
//   stx is the best step obtained so far and is an
//        endpoint of the interval that contains the minimizer.
//
//   fx is the function at stx.
//   dx is the derivative of the function at stx.
//        The derivative must be negative in the direction of
//        the step, that is, dx and stp - stx must have opposite
//        signs.
//
//   sty is the second endpoint of the interval that contains the minimizer.
//   fy is the function at sty.
//   dy is the derivative of the function at sty.
//
//   stp is the current step.
//   fp is is the function at stp
//   dp is the the derivative of the function at stp.
//
//   bracket is an bool variable which specifies if a minimizer has been
//   bracketed.
//   stpmin is a lower bound for the step.
//   stpmax is an upper bound for the step.
template <typename T>
T
step(T& stx, T& fx, T& dx, T& sty, T& fy, T& dy,
    const T stp, const T fp, const T dp, bool& bracket, const T stpmin, const T stpmax)
{
    // Some checks (they were documented).
    if (bracket)
        assert(stp > std::min(stx, sty) && stp < std::max(stx, sty));

    T stpf(0);
    const T sgnd = dp * (dx / fabs(dx));
    // First case: A higher function value. The minimum is bracketed.
    // If the cubic step is closer to stx than the quadratic step, the
    // cubic step is taken, otherwise the average of the cubic and
    // quadratic steps is taken.
    if (fp > fx) {
        const T stpc = cubic_step(stx, fx, dx, stp, fp, dp);
        const T stpq = quadratic_step(stx, fx, dx, stp, fp);
        stpf = is_closer(stpc, stpq, stx) ? stpc : (stpq + stpc) / 2;
        bracket = true;
    // Second case: A lower function value and derivatives of opposite
    // sign. The minimum is bracketed. If the cubic step is farther from
    // stp than the secant step, the cubic step is taken, otherwise the
    // secant step is taken.
    } else if (sgnd < 0) {
        const T stpc = cubic_step(stx, fx, dx, stp, fp, dp);
        const T stpq = secant_step(stx, dx, stp, dp);
        stpf = is_closer(stpc, stpq, stp) ? stpq : stpc;
        bracket = true;
    // Third case: A lower function value, derivatives of the same sign,
    // and the magnitude of the derivative decreases.
    } else if (fabs(dp) < fabs(dx)) {
        // The cubic step is computed only if the cubic tends to infinity
        // in the direction of the step or if the minimum of the cubic
        // is beyond stp. Otherwise the cubic step is defined to be the
        // secant step.
        const T theta = T(3.0) * (fx - fp) / (stp - stx) + dx + dp;
        const T s = std::max({fabs(theta), fabs(dx), fabs(dp)});
        // The case gamma = 0 only arises if the cubic does not tend
        // to infinity in the direction of the step.
        T gamma = s * sqrt(std::max(0.0, pow(theta / s, 2) - dx / s * (dp / s)));
        gamma *= signum(stx - stp);
        const T p = gamma - dp + theta;
        const T q = gamma + (dx - dp) + gamma;
        const T r = p / q;
        T stpc = stpmin;
        if (r < 0 && gamma != 0) {
            stpc = stp + r * (stx - stp);
        } else if (stp > stx) {
            stpc = stpmax;
        }
        const T stpq = secant_step(stx, dx, stp, dp);

        if (bracket) {
            // A minimizer has been bracketed. If the cubic step is
            // closer to stp than the secant step, the cubic step is
            // taken, otherwise the secant step is taken.
            stpf = is_closer(stpc, stpq, stp) ? stpc : stpq;
            if (stp > stx) {
                stpf = std::min(stp + 0.66 * (sty - stp), stpf);
            } else {
                stpf = std::max(stp + 0.66 * (sty - stp), stpf);
            }
        } else {
            // A minimizer has not been bracketed. If the cubic step is
            // farther from stp than the secant step, the cubic step is
            // taken, otherwise the secant step is taken.
            stpf = !is_closer(stpc, stpq, stp) ? stpc : stpq;
            stpf = std::min(stpmax, stpf);
            stpf = std::max(stpmin, stpf);
        }
    // Fourth case: A lower function value, derivatives of the same sign,
    // and the magnitude of the derivative does not decrease. If the
    // minimum is not bracketed, the step is either stpmin or stpmax,
    // otherwise the cubic step is taken.
    } else {
        if (bracket) {
            stpf = cubic_step(sty, fy, dy, stp, fp, dp);
        } else{
            stpf = (stp > stx) ? stpmax : stpmin;
        }
    }
    // Update the interval which contains a minimizer.
    if (fp > fx) {
        sty = stp;
        fy = fp;
        dy = dp;
    } else {
        if (sgnd < 0) {
            sty = stx;
            fy = fx;
            dy = dx;
        }
        stx = stp;
        fx = fp;
        dx = dp;
    }
    return stpf;
}

} // ns more_thuente
} // ns line_search
} // ns ook

#endif
