#ifndef LINE_SEARCH_H_
#define LINE_SEARCH_H_

#include <string>
#include <stdexcept>
#include <cassert>

enum class task_value{
    start,
    fg,
    warning_rounding_error_prevents_progress,
    warning_xtol_satisfied,
    warning_stp_eq_stpmax,
    warning_stp_eq_stpmin,
    convergence
};

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

inline
bool
sufficient_decrease_condition(const double fxap, const double fx, const double dfx_dot_p, const double a, const double c)
{
    return fxap <= fx + a * c * dfx_dot_p;
}

inline
bool
curvature_condition(const double dfxap_dot_p, const double dfx_dot_p, const double c)
{
    return fabs(dfxap_dot_p) <= c * fabs(dfx_dot_p);
}

int dcsrch(const double finit, const double ginit, double& stp, double f, double g, task_value& task, const options& opts);
double
dcstep(double& stx, double& fx, double& dx, double& sty, double& fy, double& dy, 
            const double& stp, const double& fp, const double& dp, bool& brackt, const double& stpmin, const double& stpmax);

struct dcsrch_struct{

    dcsrch_struct(double f0, double g0, double stp, double width0){
        finit = f0;
        ginit = g0;
        brackt = false;
        stage = 1;
        width = width0;
        width1 = 2.0 * width;
        /* The variables stx, fx, gx contain the values of the step,
        function, and derivative at the best step.
        The variables sty, fy, gy contain the value of the step,
        function, and derivative at sty.
        The variables stp, f, g contain the values of the step,
        function, and derivative at stp. */
        stx = 0.;
        fx = finit;
        gx = ginit;
        sty = 0.;
        fy = finit;
        gy = ginit;
        stmin = 0.;
        stmax = stp + 4.0 * stp;        
    }

    task_value
    operator()(double& stp, double f, double g, const options& opts)
    {
        validate_arguments(stp, ginit, opts);
        const bool sufficient_decrease = sufficient_decrease_condition(f, finit, ginit, stp, opts.ftol);
        const bool curvature = curvature_condition(g, ginit, opts.gtol);

        /*     Test for warnings. */
        if (brackt && (stp <= stmin || stp >= stmax)) {
            return task_value::warning_rounding_error_prevents_progress;
        }
        if (brackt && stmax - stmin <= opts.xtol * stmax) {
            return task_value::warning_xtol_satisfied;        
        }
        if (stp == opts.stpmax && sufficient_decrease && curvature) {
            return task_value::warning_stp_eq_stpmax;
        }
        if (stp == opts.stpmin && (!sufficient_decrease || !curvature)) {
            return task_value::warning_stp_eq_stpmin;
        }
        if (stage == 1 && sufficient_decrease && g >= 0.){
            stage = 2;
        }        
        /* A modified function is used to predict the step during the
        first stage if a lower function value has been obtained but
        the decrease is not sufficient. */
        if (stage == 1 && f <= fx && !sufficient_decrease) {
            /* Define the modified function and derivative values. */
            const double gtest = opts.ftol * ginit;                
            const double fm = f - stp * gtest;
            double fxm = fx - stx * gtest;
            double fym = fy - sty * gtest;
            const double gm = g - gtest;
            double gxm = gx - gtest;
            double gym = gy - gtest;
            /* Call dcstep to update stx, sty, and to compute the new step. */
            stp = dcstep(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, brackt, stmin, stmax);
            /* Reset the function and derivative values for f. */
            fx = fxm + stx * gtest;
            fy = fym + sty * gtest;
            gx = gxm + gtest;
            gy = gym + gtest;
        } else  {
            /* Call dcstep to update stx, sty, and to compute the new step. */
            stp = dcstep(stx, fx, gx, sty, fy, gy, stp, f, g, brackt, stmin, stmax);
        }
        /* Decide if a bisection step is needed. */
        if (brackt) {
            if (fabs(sty - stx) >= 0.66 * width1) {
                stp = stx + 0.5 * (sty - stx);
            }
            width1 = width;
            width = fabs(sty - stx);

            /* Set the minimum and maximum steps allowed for stp. */
            stmin = std::min(stx, sty);
            stmax = std::max(stx, sty);
        } else {
            stmin = stp + 1.1 * (stp - stx);
            stmax = stp + 4.0 * (stp - stx);
        }
        /* Force the step to be within the bounds opts.stpmax and opts.stpmin. */
        stp = std::max(stp, opts.stpmin);
        stp = std::min(stp, opts.stpmax);
        /* If further progress is not possible, let stp be the best point obtained
        during the search. */
        if ((brackt && (stp <= stmin || stp >= stmax) )
            || (brackt && stmax - stmin <= opts.xtol * stmax))
        {
            stp = stx;
        }
        /* Obtain another function and derivative. */
        return task_value::fg;        
    }

    int stage;
    double finit;
    double ginit;
    double width;
    double stmin;
    double stmax;
    double width1;
    double fx;
    double fy;
    double gx;
    double gy;
    bool brackt;
    double stx, sty;
};

#endif
