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

#ifndef OOK_LINE_SEARCH_MORE_THUENTE_MCSTEP_HPP_
#define OOK_LINE_SEARCH_MORE_THUENTE_MCSTEP_HPP_

#include "ook/line_search/step_functions.hpp"

namespace ook{
namespace line_search{

/// \brief The purpose of mcstep is to compute a safeguarded step for a
/// linesearch and to update an interval of uncertainty for a minimizer of the
/// function.
/// \detail The parameter stx contains the step with the least function value.
/// The parameter stp contains the current step. It is assumed that the
/// derivative at stx is negative in the direction of the step. If brackt is
/// set to true, then a minimizer has been bracketed in an interval of
/// uncertainty with endpoints stx, sty.
///
/// \param stx The step with the least function value.
/// \param fx The best function value obtained so far.
/// \param dx The derivative at the best point obtained so far.
/// \param sty The step at the endpoint of the interval.
/// \param fy The function value at the endpoint of the interval.
/// \param dy The derivative at the endpoint of the interval.
/// \param stp The current step.
/// \param fp The function value at the current step.
/// \param dp The derivative at the current step.
/// \param brackt Boolean specifying if a minimizer has been bracketed.
///        If brackt is true, then on input stp must be between stx and sty.
///        On output, stp is set to the new step.
/// \param stpmin Lower bound for the step.
/// \param stpmax Upper bound for the step.
/// \return An integer in [0,1,..,5] indicating which step has been used, or
/// zero if there is a problem with the input parameters.
template <typename T>
int mcstep(T& stx, T& fx, T& dx, T& sty, T& fy,
            T& dy, T& stp, T& fp, T& dp, bool& brackt, T& stpmin, T& stpmax)
{
    int info = 0;
    bool invalid_bracket = (stp <= std::min(stx, sty)
                            || stp >= std::max(stx, sty));

    if((brackt && invalid_bracket) || dx * (stp - stx) >= 0.0)
    {
        return info;
    }

    T stpf(0);
    bool bound;

    // Determine if the derivatives have opposite sign.
    const T sgnd = dp * (dx / fabs(dx));

    // First case: A higher function value. The minimum is bracketed. If the
    // cubic step is closer to stx than the quadratic step, the cubic step is
    // taken, otherwise the average of the cubic and quadratic steps is taken.
    if(fp > fx){
        info = 1;
        bound = true;
        const T stpc = cubic_step(stx, fx, dx, stp, fp, dp);
        const T stpq = quadratic_step(stx, fx, dx, stp, fp);
        stpf = is_closer(stpc, stpq, stx) ? stpc : (stpq + stpc) / 2;
        brackt = true;
    // Second case: A lower function value and derivatives of opposite sign.
    // The minimum is bracketed. If the cubic step is farther from stp than
    // the secant step, the cubic step is taken, otherwise the secant step is
    // taken.
    } else if(sgnd < 0.0){
        info = 2;
        bound = false;
        const T stpc = cubic_step(stx, fx, dx, stp, fp, dp);
        const T stpq = secant_step(stx, dx, stp, dp);
        stpf = is_closer(stpc, stpq, stp) ? stpq : stpc;
        brackt = true;
    // Third case: A lower function value, derivatives of the same sign,and the
    // magnitude of the derivative decreases. The cubic step is computed only
    // if the cubic tends to infinity in the direction of the step or if the
    // minimum of the cubic is beyond stp. Otherwise the cubic step is defined
    // to be the secant step.
    } else if(fabs(dp) < fabs(dx)){
        info = 3;
        bound = true;
        const T theta = 3.0 * (fx - fp)  / (stp - stx) + dx + dp;
        const T s = std::max({fabs(theta), fabs(dx), fabs(dp)});
        // The case gamma = 0 only arises if the cubic does not tend to
        // infinity in the direction of the step.
        const T ts = theta / s;
        T gamma = s * sqrt(std::max(0.0, ts * ts - dx / s * (dp / s)));
        if(stp > stx){
            gamma = -gamma;
        }
        const T p = gamma - dp + theta;
        const T q = gamma + (dx - dp) + gamma;
        const T r = p / q;

        T stpc = stpmin;
        if(r < 0.0 && gamma != 0.0){
            stpc = stp + r * (stx - stp);
        } else if(stp > stx){
            stpc = stpmax;
        }

        const T stpq = stp + dp / (dp - dx) * (stx - stp);
        if(brackt){
            // A minimizer has been bracketed. If the cubic step is closer to
            // stp than the secant step, the cubic step is taken, otherwise the
            // secant step is taken.
            stpf = is_closer(stpc, stpq, stp) ? stpc : stpq;
        }else{
            // A minimizer has not been bracketed. If the cubic step is farther
            // from stp than the secant step, the cubic step is taken,
            // otherwise the secant step is taken.
            stpf = !is_closer(stpc, stpq, stp) ? stpc : stpq;
        }
    // Fourth case: A lower function value, derivatives of the same sign, and
    // the magnitude of the derivative does not decrease. If the minimum is not
    // bracketed, the step is either stpmin or stpmax, otherwise the cubic step
    //  is taken.
    } else {
        info = 4;
        bound = false;
        if(brackt){
            stpf = cubic_step(sty, fy, dy, stp, fp, dp);
        }else{
            stpf = (stp > stx) ? stpmax : stpmin;
        }
    }

    // Update the interval of uncertainty. This update does not depend on the
    // new step or the case analysis above.
    if(fp > fx) {
        sty = stp;
        fy = fp;
        dy = dp;
    } else {
        if(sgnd < 0.f) {
            sty = stx;
            fy = fx;
            dy = dx;
        }
        stx = stp;
        fx = fp;
        dx = dp;
    }
    // Compute the new step and safeguard it
    stpf = std::min(stpmax, stpf);
    stpf = std::max(stpmin, stpf);
    stp = stpf;
    if(brackt && bound){
        if(sty > stx){
            stp = std::min(stx + (sty - stx) * 0.66, stp);
        }else{
            stp = std::max(stx + (sty - stx) * 0.66, stp);
        }
    }
    return info;
}

} // ns line_search
} // ns ook

#endif
