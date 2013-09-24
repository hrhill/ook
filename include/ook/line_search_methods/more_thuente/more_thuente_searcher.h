#ifndef OOK_LINE_SEARCH_MORE_THUENTE_SEARCHER_H_
#define OOK_LINE_SEARCH_MORE_THUENTE_SEARCHER_H_

#include <cmath>

#include "../../state_value.h"

#include "line_search_conditions.h"
#include "safe_step.h"

namespace ook{

namespace line_search{

template <typename T>
struct more_thuente_searcher{

    more_thuente_searcher(T f0_, T g0_, T stp, T width0)
    :
        f0(f0_), g0(g0_), brackt(false), stage(1), width(width0),
            width1(2.0 * width0), 
                stx(0), fx(f0), gx(g0), 
                    sty(0), fy(f0), gy(0),
                        stmin(0), stmax(5.0 * stp)
    {}

    template <typename Options>
    std::pair<state_value, T>
    operator()(T stp, T f, T g, const Options& opts)
    {
        const bool sufficient_decrease = sufficient_decrease_condition(f0, f, g0, stp, opts.ftol);
        const bool curvature = curvature_condition(g0, g, opts.gtol);

        /*     Test for warnings. */
        if (brackt && (stp <= stmin || stp >= stmax)) {
            return std::make_pair(state_value::warning_rounding_error_prevents_progress, stp);
        }
        if (brackt && stmax - stmin <= opts.xtol * stmax) {
            return std::make_pair(state_value::warning_xtol_satisfied, stp);        
        }
        if (stp == opts.stpmax && sufficient_decrease && curvature) {
            return std::make_pair(state_value::warning_stp_eq_stpmax, stp);
        }
        if (stp == opts.stpmin && (!sufficient_decrease || !curvature)) {
            return std::make_pair(state_value::warning_stp_eq_stpmin, stp);
        }
        if (stage == 1 && sufficient_decrease && g >= T(0.0)){
            stage = 2;
        }        
        /* A modified function is used to predict the step during the
        first stage if a lower function value has been obtained but
        the decrease is not sufficient. */
        if (stage == 1 && f <= fx && !sufficient_decrease) {
            // Define the modified function and derivative values. 
            const T gtest = opts.ftol * g0;                
            const T fm = f - stp * gtest;
            T fxm = fx - stx * gtest;
            T fym = fy - sty * gtest;
            const T gm = g - gtest;
            T gxm = gx - gtest;
            T gym = gy - gtest;
            // Call safe_step to update stx, sty, and to compute the new step. 
            stp = safe_step(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, brackt, stmin, stmax);
            // Reset the function and derivative values for f.
            fx = fxm + stx * gtest;
            fy = fym + sty * gtest;
            gx = gxm + gtest;
            gy = gym + gtest;
        } else  {
            // Call safe_step to update stx, sty, and to compute the new step.
            stp = safe_step(stx, fx, gx, sty, fy, gy, stp, f, g, brackt, stmin, stmax);
        }
        /* Decide if a bisection step is needed. */
        if (brackt) {
            if (fabs(sty - stx) >= T(0.66) * width1) {
                stp = stx + T(0.5) * (sty - stx);
            }
            width1 = width;
            width = fabs(sty - stx);

            // Set the minimum and maximum steps allowed for stp.
            stmin = std::min(stx, sty);
            stmax = std::max(stx, sty);
        } else {
            stmin = stp + T(1.1) * (stp - stx);
            stmax = stp + T(4.0) * (stp - stx);
        }
        // Force the step to be within the bounds opts.stpmax and opts.stpmin.
        stp = std::max(stp, opts.stpmin);
        stp = std::min(stp, opts.stpmax);
        //If further progress is not possible, let stp be the best point obtained during the search.
        if ((brackt && (stp <= stmin || stp >= stmax))
            || (brackt && stmax - stmin <= opts.xtol * stmax))
        {
            stp = stx;
        }
        return std::make_pair(state_value::update, stp);        
    }

    const T f0;
    const T g0;
    bool brackt;    
    int stage;    
    T width;
    T width1;
    T stx;    
    T fx;
    T gx;
    T sty;    
    T fy;
    T gy;
    T stmin;
    T stmax;
};

} // ns line search

} // ns ook

#endif