#include <stdbool.h>
#include <string.h>
#include <math.h>

#include "line_search.h"

#define max(x, y) ((x) > (y) ? (x) : (y))
#define min(x, y) ((x) < (y) ? (x) : (y))

typedef short ftnlen;

int s_cmp(char *a0, char *b0, ftnlen la, ftnlen lb)
{
    unsigned char *a, *aend, *b, *bend;
    a = (unsigned char *)a0;
    b = (unsigned char *)b0;
    aend = a + la;
    bend = b + lb;

    if(la <= lb){
        while(a < aend)
            if(*a != *b)
                return( *a - *b );
            else{ 
                ++a; ++b; 
            }

        while(b < bend)
            if(*b != ' ')
                return( ' ' - *b );
            else    ++b;
    }else{
        while(b < bend)
            if(*a == *b)
                { ++a; ++b; }
            else
                return( *a - *b );
        while(a < aend)
            if(*a != ' ')
                return(*a - ' ');
            else    ++a;
    }
    return(0);
}

int dcsrch_(double* stp, double* f, double* g, double* ftol, double* gtol, double* xtol, std::string& task,  double* stpmin, double* stpmax, ftnlen task_len)
{
/*  Subroutine dcsrch

    This subroutine finds a step that satisfies a sufficient
    decrease condition and a curvature condition.

    Each call of the subroutine updates an interval with
    endpoints stx and sty. The interval is initially chosen
    so that it contains a minimizer of the modified function

        psi(stp) = f(stp) - f(0) - ftol*stp*f'(0).

    If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the
    interval is chosen so that it contains a minimizer of f.

    The algorithm is designed to find a step that satisfies
    the sufficient decrease condition

        f(stp) <= f(0) + ftol*stp*f'(0),

    and the curvature condition

        abs(f'(stp)) <= gtol*abs(f'(0)).

    If ftol is less than gtol and if, for example, the function
    is bounded below, then there is always a step which satisfies
    both conditions.

    If no step can be found that satisfies both conditions, then
    the algorithm stops with a warning. In this case stp only
    satisfies the sufficient decrease condition.

    A typical invocation of dcsrch has the following outline:

    Evaluate the function at stp = 0.0d0; store in f.
    Evaluate the gradient at stp = 0.0d0; store in g.
    Choose a starting step stp.

    task = 'START'
    10 continue
        call dcsrch(stp,f,g,ftol,gtol,xtol,task = stpmin,stpmax,isave,dsave)
        if (task .eq. 'FG') then
            Evaluate the function and the gradient at stp
            go to 10
        end if

    NOTE: The user must not alter work arrays between calls.

    The subroutine statement is

    subroutine dcsrch(f,g,stp,ftol,gtol,xtol,stpmin,stpmax,task = isave,dsave)
    where

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

        ftol is a double precision variable.
            On entry ftol specifies a nonnegative tolerance for the sufficient decrease condition.
            On exit ftol is unchanged.

        gtol is a double precision variable.
            On entry gtol specifies a nonnegative tolerance for the curvature condition.
            On exit gtol is unchanged.

        xtol is a double precision variable.
            On entry xtol specifies a nonnegative relative tolerance
                for an acceptable step. The subroutine exits with a
                warning if the relative difference between sty and stx
                is less than xtol.
            On exit xtol is unchanged.

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

        stpmin is a double precision variable.
            On entry stpmin is a nonnegative lower bound for the step.
            On exit stpmin is unchanged. 

        stpmax is a double precision variable.
            On entry stpmax is a nonnegative upper bound for the step. 
            On exit stpmax is unchanged.

        MINPACK-1 Project. June 1983.
        Argonne National Laboratory.
        Jorge J. More' and David J. Thuente.

        MINPACK-2 Project. November 1993.
        Argonne National Laboratory and University of Minnesota. 
        Brett M. Averick, Richard G. Carter, and Jorge J. More'. */

    static int stage;
    static double finit, ginit, width, ftest, gtest, stmin, stmax, width1, fm, gm, fx, fy, gx, gy;
    static bool brackt;
    static double fxm, fym, gxm, gym, stx, sty;

    /* Function Body */
    if (task == "START") {
        /*        Check the input arguments for errors. */
        if (*stp < *stpmin) {
            task = "ERROR: STP .LT. STPMIN";
        }
        if (*stp > *stpmax) {
            task = "ERROR: STP .GT. STPMAX";
        }
        if (*g >= 0.) {
            task = "ERROR: INITIAL G .GE. ZERO";
        }
        if (*ftol < 0.) {
            task = "ERROR: FTOL .LT. ZERO";
        }
        if (*gtol < 0.) {
            task = "ERROR: GTOL .LT. ZERO";
        }
        if (*xtol < 0.) {
            task = "ERROR: XTOL .LT. ZERO";
        }
        if (*stpmin < 0.) {
            task = "ERROR: STPMIN .LT. ZERO";
        }
        if (*stpmax < *stpmin) {
            task = "ERROR: STPMAX .LT. STPMIN";
        }
        /* Exit if there are errors on input. */
        if (task.find("ERROR") != std::string::npos) {
            return 0;
        }
        /*        Initialize local variables. */
        brackt = false;
        stage = 1;
        finit = *f;
        ginit = *g;
        gtest = *ftol * ginit;
        width = *stpmax - *stpmin;
        width1 = width / .5;
        /*        The variables stx, fx, gx contain the values of the step, */
        /*        function, and derivative at the best step. */
        /*        The variables sty, fy, gy contain the value of the step, */
        /*        function, and derivative at sty. */
        /*        The variables stp, f, g contain the values of the step, */
        /*        function, and derivative at stp. */
        stx = 0.;
        fx = finit;
        gx = ginit;
        sty = 0.;
        fy = finit;
        gy = ginit;
        stmin = 0.;
        stmax = *stp + *stp * 4.;

        task =  "FG";
        return 0;
    }
    /*     If psi(stp) <= 0 and f'(stp) >= 0 for some step, then the */
    /*     algorithm enters the second stage. */
    ftest = finit + *stp * gtest;
    if (stage == 1 && *f <= ftest && *g >= 0.) {
        stage = 2;
    }
    /*     Test for warnings. */
    if (brackt && (*stp <= stmin || *stp >= stmax)) {
        task = "WARNING: ROUNDING ERRORS PREVENT PROGRESS";
    }
    if (brackt && stmax - stmin <= *xtol * stmax) {
        task = "WARNING: XTOL TEST SATISFIED";
    }
    if (*stp == *stpmax && *f <= ftest && *g <= gtest) {
        task = "WARNING: STP = STPMAX";
    }
    if (*stp == *stpmin && (*f > ftest || *g >= gtest)) {
        task = "WARNING: STP = STPMIN";
    }
    /*     Test for convergence. */
    if (*f <= ftest && fabs(*g) <= *gtol * (-ginit)) {
        task = "CONVERGENCE";
    }
    /*     Test for termination. */
    if (task.find("WARN") != std::string::npos || task.find("CONV") != std::string::npos) {
        return 0;
    }
    /*     A modified function is used to predict the step during the */
    /*     first stage if a lower function value has been obtained but */
    /*     the decrease is not sufficient. */
    if (stage == 1 && *f <= fx && *f > ftest) {
    /*        Define the modified function and derivative values. */
        fm = *f - *stp * gtest;
        fxm = fx - stx * gtest;
        fym = fy - sty * gtest;
        gm = *g - gtest;
        gxm = gx - gtest;
        gym = gy - gtest;
        /*        Call dcstep to update stx, sty, and to compute the new step. */
        dcstep_(&stx, &fxm, &gxm, &sty, &fym, &gym, stp, &fm, &gm, &brackt, &stmin, &stmax);
        /*        Reset the function and derivative values for f. */
        fx = fxm + stx * gtest;
        fy = fym + sty * gtest;
        gx = gxm + gtest;
        gy = gym + gtest;
    } else {
        /*       Call dcstep to update stx, sty, and to compute the new step. */
        dcstep_(&stx, &fx, &gx, &sty, &fy, &gy, stp, f, g, &brackt, &stmin, &stmax);
    }
    /*     Decide if a bisection step is needed. */
    if (brackt) {
        if (fabs(sty - stx) >= width1 * .66) {
            *stp = stx + (sty - stx) * .5;
        }
        width1 = width;
        width = fabs(sty - stx);
    }
    /*     Set the minimum and maximum steps allowed for stp. */
    if (brackt) {
        stmin = min(stx,sty);
        stmax = max(stx,sty);
    } else {
        stmin = *stp + (*stp - stx) * 1.1;
        stmax = *stp + (*stp - stx) * 4.;
    }
    /*     Force the step to be within the bounds stpmax and stpmin. */
    *stp = max(*stp, *stpmin);
    *stp = min(*stp, *stpmax);
    /*     If further progress is not possible, let stp be the best */
    /*     point obtained during the search. */
    if (brackt && (*stp <= stmin || *stp >= stmax) || brackt && stmax - stmin <= *xtol * stmax){
        *stp = stx;
    }
    /*     Obtain another function and derivative. */
    task = "FG";
} 
