#include <cmath>
#include <algorithm>
#include <stdexcept>

#include "line_search.h"

#define ASSERT_CHECK_AND_THROW(predicate) assert((predicate));\
    if (!(predicate))\
        throw std::invalid_argument(std::string(#predicate) + std::string(" must hold."));

void
validate_arguments(const double stp, const double g, const options& opts)
{
    ASSERT_CHECK_AND_THROW(!(stp < opts.stpmin));
    ASSERT_CHECK_AND_THROW(!(stp > opts.stpmax));
    ASSERT_CHECK_AND_THROW(!(g > 0.0));
}

#undef ASSERT_CHECK_AND_THROW

int dcsrch(const double finit, const double ginit, double& stp, double f, double g, task_value& task, const options& opts)
{
/*  Subroutine dcsrch

    This subroutine finds a step that satisfies a sufficient
    decrease condition and a curvature condition.

    Each call of the subroutine updates an interval with
    endpoints stx and sty. The interval is initially chosen
    so that it contains a minimizer of the modified function

        psi(stp) = f(stp) - f(0) - opts.ftolstpf'(0).

    If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
    interval is chosen so that it contains a minimizer of f.

    The algorithm is designed to find a step that satisfies
    the sufficient decrease condition

        f(stp) <= f(0) + opts.ftolstpf'(0),

    and the curvature condition

        abs(f'(stp)) <= opts.gtol*abs(f'(0)).

    If opts.ftol is less than opts.gtol and if, for example, the function
    is bounded below, then there is always a step which satisfies
    both conditions.

    If no step can be found that satisfies both conditions, then
    the algorithm stops with a warning. In this case stp only
    satisfies the sufficient decrease condition.

    A typical invocation of dcsrch has the following outline:

    Evaluate the function at stp = 0.0d0; store in f.
    Evaluate the gradient at stp = 0.0d0; store in g.
    Choose a starting step stp.

    task = "START"
    do{
        dcsrch(stp, f, g, task, opts);
        if (task.contains("FG")n
            Evaluate the function and the gradient at stp
    }

    Arguments:
        stp is a double precision variable.
            On entry stp is the current estimate of a satisfactory
                step. On initial entry, a positive initial estimate
                must be provided.
            On exit stp is the current estimate of a satisfactory step
                if task = 'FG'. If task = 'CONV' then stp satisfies
                the sufficient decrease and curvature condition.

        f is a double precision variable.
            On initial entry f is the value of the function at 0.
                On subsequent entries f is the value of the function at stp.
            On exit f is the value of the function at stp.

        g is a double precision variable.
            On initial entry g is the derivative of the function at 0.
            On subsequent entries g is the derivative of the function at stp.
            On exit g is the derivative of the function at stp.

        task is a character variable of length at least 60.
            On initial entry task must be set to 'START'.
            On exit task indicates the required action:

            If task(1:2) = 'FG' then evaluate the function and
                derivative at stp and call dcsrch again.

            If task(1:4) = 'CONV' then the search is successful.

            If task(1:4) = 'WARN' then the subroutine is not able
                to satisfy the convergence conditions. The exit value of
                stp contains the best point found during the search.

            If task(1:5) = 'ERROR' then there is an error in the
                input arguments.

            On exit with convergence, a warning or an error, the
                variable task contains additional information.

        opts.ftol specifies a nonnegative tolerance for the sufficient decrease condition.

        opts.gtol specifies a nonnegative tolerance for the curvature condition.

        opts.xtol specifies a nonnegative relative tolerance for an acceptable step. The 
                subroutine exits with a warning if the relative difference between sty and stx
                is less than opts.xtol.

        opts.stpmin is a nonnegative lower bound for the step.

        opts.stpmax is a nonnegative upper bound for the step. 
*/

    static int stage;
    static double width, stmin, stmax, width1, fm, gm, fx, fy, gx, gy;
    static bool brackt;
    static double stx, sty;

    validate_arguments(stp, ginit, opts);
    /* Function Body */
    if (task == task_value::start) {
        /*  Initialize local variables. */
        brackt = false;
        stage = 1;
        width = opts.stpmax - opts.stpmin;
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

        task =  task_value::fg;
        return 0;
    }
    /* If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
    algorithm enters the second stage. */
    const double gtest = opts.ftol * ginit;        
    const double ftest = finit + stp * gtest;

    if (stage == 1 && f <= ftest && g >= 0.){
        stage = 2;
    }
    /*     Test for warnings. */
    if (brackt && (stp <= stmin || stp >= stmax)) {
        task = task_value::warning_rounding_error_prevents_progress;
        return 0;        
    }
    if (brackt && stmax - stmin <= opts.xtol * stmax) {
        task = task_value::warning_xtol_satisfied;        
        return 0;        
    }
    if (stp == opts.stpmax && f <= ftest && g <= gtest) {
        task = task_value::warning_stp_eq_stpmax;
        return 0;        
    }
    if (stp == opts.stpmin && (f > ftest || g >= gtest)) {
        task = task_value::warning_stp_eq_stpmin;
        return 0;
    }
    /*     Test for convergence. */
    if (f <= ftest && fabs(g) <= opts.gtol * (-ginit)) {
        task = task_value::convergence;
        return 0;
    }
    /* A modified function is used to predict the step during the
    first stage if a lower function value has been obtained but
    the decrease is not sufficient. */
    if (stage == 1 && f <= fx && f > ftest) {
        /* Define the modified function and derivative values. */
        double fm = f - stp * gtest;
        double fxm = fx - stx * gtest;
        double fym = fy - sty * gtest;
        double gm = g - gtest;
        double gxm = gx - gtest;
        double gym = gy - gtest;
        /* Call dcstep to update stx, sty, and to compute the new step. */
        dcstep(stx, fxm, gxm, sty, fym, gym, stp, fm, gm, brackt, stmin, stmax);
        /* Reset the function and derivative values for f. */
        fx = fxm + stx * gtest;
        fy = fym + sty * gtest;
        gx = gxm + gtest;
        gy = gym + gtest;
    } else {
        /* Call dcstep to update stx, sty, and to compute the new step. */
        dcstep(stx, fx, gx, sty, fy, gy, stp, f, g, brackt, stmin, stmax);
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
    if (brackt && (stp <= stmin || stp >= stmax) 
        || brackt && stmax - stmin <= opts.xtol * stmax)
    {
        stp = stx;
    }
    /*     Obtain another function and derivative. */
    task = task_value::fg;
} 
