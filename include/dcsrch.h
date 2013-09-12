#ifndef OOK_DCSRCH_H_
#define OOK_DCSRCH_H_

#include <cmath>
#include "line_search_conditions.h"
#include "dcstep.h"
#include "state.h"

namespace ook{

#define ASSERT_CHECK_AND_THROW(predicate) assert((predicate));\
    if (!(predicate))\
        throw std::invalid_argument(std::string(#predicate) + std::string(" must hold."));

struct options{
    options(const double ftolIn, const double gtolIn, const double xtolIn, const double stpminIn, const double stpmaxIn)
    :
        ftol(ftolIn), 
        gtol(gtolIn), 
        xtol(xtolIn), 
        stpmin(stpminIn), 
        stpmax(stpmaxIn)
    {
        /* Check the input arguments for errors. */
        ASSERT_CHECK_AND_THROW(ftol >= 0.);
        ASSERT_CHECK_AND_THROW(gtol >= 0.); 
        ASSERT_CHECK_AND_THROW(xtol >= 0.);
        ASSERT_CHECK_AND_THROW(stpmin >= 0.); 
        ASSERT_CHECK_AND_THROW(stpmax > stpmin); 
    }
    const double ftol;
    const double gtol;
    const double xtol;
    const double stpmin;
    const double stpmax;
};

inline
void
validate_arguments(const double stp, const double g, const options& opts)
{
    ASSERT_CHECK_AND_THROW(!(stp < opts.stpmin));
    ASSERT_CHECK_AND_THROW(!(stp > opts.stpmax));
    ASSERT_CHECK_AND_THROW(!(g > 0.0));
}

#undef ASSERT_CHECK_AND_THROW

struct dcsrch_struct{

    dcsrch_struct(double f0_, double g0_, double stp, double width0)
    :
        f0(f0_), g0(g0_), brackt(false), stage(1), width(width0),
            width1(2.0 * width0), 
                stx(0), fx(f0), gx(g0), 
                    sty(0), fy(f0), gy(0),
                        stmin(0), stmax(5.0 * stp)
    {}

    std::pair<state_value, double>
    operator()(double stp, double f, double g, const options& opts)
    {
        validate_arguments(stp, g0, opts);
        const bool sufficient_decrease = sufficient_decrease_condition(state(f0, f, g0, g, stp, state_value::update), opts.ftol);
        const bool curvature = auxilliary_curvature_condition(state(f0, f, g0, g, stp, state_value::update), opts.gtol);

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
        if (stage == 1 && sufficient_decrease && g >= 0.){
            stage = 2;
        }        
        /* A modified function is used to predict the step during the
        first stage if a lower function value has been obtained but
        the decrease is not sufficient. */
        if (stage == 1 && f <= fx && !sufficient_decrease) {
            // Define the modified function and derivative values. 
            const double gtest = opts.ftol * g0;                
            const double fm = f - stp * gtest;
            double fxm = fx - stx * gtest;
            double fym = fy - sty * gtest;
            const double gm = g - gtest;
            double gxm = gx - gtest;
            double gym = gy - gtest;
            // Call dcstep to update stx, sty, and to compute the new step. 
            stp = dcstep(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, brackt, stmin, stmax);
            // Reset the function and derivative values for f.
            fx = fxm + stx * gtest;
            fy = fym + sty * gtest;
            gx = gxm + gtest;
            gy = gym + gtest;
        } else  {
            // Call dcstep to update stx, sty, and to compute the new step.
            stp = dcstep(stx, fx, gx, sty, fy, gy, stp, f, g, brackt, stmin, stmax);
        }
        /* Decide if a bisection step is needed. */
        if (brackt) {
            if (fabs(sty - stx) >= 0.66 * width1) {
                stp = stx + 0.5 * (sty - stx);
            }
            width1 = width;
            width = fabs(sty - stx);

            // Set the minimum and maximum steps allowed for stp.
            stmin = std::min(stx, sty);
            stmax = std::max(stx, sty);
        } else {
            stmin = stp + 1.1 * (stp - stx);
            stmax = stp + 4.0 * (stp - stx);
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

    const double f0;
    const double g0;
    bool brackt;    
    int stage;    
    double width;
    double width1;
    double stx;    
    double fx;
    double gx;
    double sty;    
    double fy;
    double gy;
    double stmin;
    double stmax;

};

}

#endif