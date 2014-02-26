#ifndef MCSRCH_H_
#define MCSRCH_H_

#include "../line_search/more_thuente/safe_step.h"

struct lb3_1_ {
    int mp, lp;
    double gtol, stpmin, stpmax;
} lb3_1 = { 6, 6, .9, 1e-20, 1e20 };


template <typename T>
int mcsrch(int n, T *x, T f, T dg, T *s, T& stp, T ftol, T xtol, int maxfev, int& info, int& nfev, T *wa)
{
/*                     SUBROUTINE MCSRCH */

/*     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES */
/*     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION. */

/*     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF */
/*     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF */
/*     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A */
/*     MINIMIZER OF THE MODIFIED FUNCTION */

/*          F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S). */

/*     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION */
/*     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE, */
/*     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT */
/*     CONTAINS A MINIMIZER OF F(X+STP*S). */

/*     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES */
/*     THE SUFFICIENT DECREASE CONDITION */

/*           F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S), */

/*     AND THE CURVATURE CONDITION */

/*           fabs(GRADF(X+STP*S)'S)) .LE. GTOL*fabs(GRADF(X)'S). */

/*     IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION */
/*     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES */
/*     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH */
/*     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING */
/*     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY */
/*     SATISFIES THE SUFFICIENT DECREASE CONDITION. */

/*     THE SUBROUTINE STATEMENT IS */

/*        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA) */
/*     WHERE */

/*       N IS A POSITIVE int INPUT VARIABLE SET TO THE NUMBER */
/*         OF VARIABLES. */

/*       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE */
/*         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS */
/*         X + STP*S. */

/*       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F */
/*         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S. */

/*       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE */
/*         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT */
/*         OF F AT X + STP*S. */

/*       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE */
/*         SEARCH DIRECTION. */

/*       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN */
/*         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT */
/*         STP CONTAINS THE FINAL ESTIMATE. */

/*       FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse */
/*         communication implementation GTOL is defined in a COMMON */
/*         statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE */
/*         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE */
/*         SATISFIED. */

/*       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS */
/*         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY */
/*         IS AT MOST XTOL. */

/*       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH */
/*         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse */
/*         communication implementatin they are defined in a COMMON */
/*         statement). */

/*       MAXFEV IS A POSITIVE int INPUT VARIABLE. TERMINATION */
/*         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST */
/*         MAXFEV BY THE END OF AN ITERATION. */

/*       INFO IS AN int OUTPUT VARIABLE SET AS FOLLOWS: */

/*         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT. */

/*         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE */
/*                   DIRECTIONAL DERIVATIVE CONDITION HOLD. */

/*         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY */
/*                   IS AT MOST XTOL. */

/*         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV. */

/*         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN. */

/*         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX. */

/*         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS. */
/*                   THERE MAY NOT BE A STEP WHICH SATISFIES THE */
/*                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS. */
/*                   TOLERANCES MAY BE TOO SMALL. */

/*       NFEV IS AN int OUTPUT VARIABLE SET TO THE NUMBER OF CALLS TO FCN. */
/*       WA IS A WORK ARRAY OF LENGTH N. */

    const double p5 = .5;
    const double p66 = .66;
    const double xtrapf = 4.;

    /* Local variables */
    static double finit, width, stmin, stmax;
    static bool stage1;
    static double width1, ftest1, fx, fy;
    static bool brackt;
    static double dginit, dgtest;
    static double dgx, dgy, stx, sty;
    if (info == 0) {
        /*     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION AND CHECK THAT S IS A DESCENT DIRECTION. */
        dginit = dg;

        /*  INITIALIZE LOCAL VARIABLES. */
        brackt = false;
        stage1 = true;
        nfev = 0;
        finit = f;
        dgtest = ftol * dginit;
        width = lb3_1.stpmax - lb3_1.stpmin;
        width1 = width / p5;
        std::copy(x, x + n, wa);

        /*     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP, */
        /*     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP. */
        /*     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP, */
        /*     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF */
        /*     THE INTERVAL OF UNCERTAINTY. */
        /*     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP, */
        /*     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP. */
        stx = 0.0;
        fx = finit;
        dgx = dginit;
        sty = 0.0;
        fy = finit;
        dgy = dginit;
    }

    while(true){
        if (info == 0){
            /*        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND */
            /*        TO THE PRESENT INTERVAL OF UNCERTAINTY. */
            if (brackt) {
                stmin = std::min(stx, sty);
                stmax = std::max(stx, sty);
            } else {
                stmin = stx;
                stmax = stp + xtrapf * (stp - stx);
            }
            /*  FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN. */
            stp = std::max(stp, lb3_1.stpmin);
            stp = std::min(stp, lb3_1.stpmax);

            /* IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET */
            /* STP BE THE LOWEST POINT OBTAINED SO FAR. */
            if ((brackt && (stp <= stmin || stp >= stmax)) || nfev >= maxfev - 1 || (brackt && stmax - stmin <= xtol * stmax)) {
                stp = stx;
            }
            /*        EVALUATE THE FUNCTION AND GRADIENT AT STP */
            /*        AND COMPUTE THE DIRECTIONAL DERIVATIVE. */
            for (int j = 0; j < n; ++j) {
                x[j] = wa[j] + stp * s[j];
            }
            info = -1;
            return 0;
        }
        info = 0;
        ++nfev;
        ftest1 = finit + stp * dgtest;

        if ((brackt && (stp <= stmin || stp >= stmax))) {
            info = 6;
        }
        if (stp == lb3_1.stpmax && f <= ftest1 && dg <= dgtest) {
            info = 5;
        }
        if (stp == lb3_1.stpmin && (f > ftest1 || dg >= dgtest)) {
            info = 4;
        }
        if (nfev >= maxfev) {
            info = 3;
        }
        if (brackt && stmax - stmin <= xtol * stmax) {
            info = 2;
        }
        if (f <= ftest1 && fabs(dg) <= lb3_1.gtol * (-dginit)) {
            info = 1;
        }
        if (info != 0) {
            return 0;
        }
        /*        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED */
        /*        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE. */
        if (stage1 && f <= ftest1 && dg >= std::min(ftol, lb3_1.gtol) * dginit) {
            stage1 = false;
        }
        /*        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF */
        /*        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED */
        /*        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE */
        /*        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN */
        /*        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT. */
        if (stage1 && f <= fx && f > ftest1) {
            /*           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES. */
            double fm = f - stp * dgtest;
            double fxm = fx - stx * dgtest;
            double fym = fy - sty * dgtest;
            double dgm = dg - dgtest;
            double dgxm = dgx - dgtest;
            double dgym = dgy - dgtest;
            /* CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY AND TO COMPUTE THE NEW STEP. */
            //stp = mcstep(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt, stmin, stmax);
            stp = ook::line_search::safe_step(stx, fxm, dgxm, sty, fym, dgym, stp, fm, dgm, brackt, stmin, stmax);
            /* RESET THE FUNCTION AND GRADIENT VALUES FOR F. */
            fx = fxm + stx * dgtest;
            fy = fym + sty * dgtest;
            dgx = dgxm + dgtest;
            dgy = dgym + dgtest;
        } else {
            /* CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY AND TO COMPUTE THE NEW STEP. */
            //stp = mcstep(stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin, stmax);
            stp = ook::line_search::safe_step(stx, fx, dgx, sty, fy, dgy, stp, f, dg, brackt, stmin, stmax);
        }
        /* FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE INTERVAL OF UNCERTAINTY. */
        if (brackt) {
            if (fabs(sty - stx) >= p66 * width1) {
                stp = stx + p5 * (sty - stx);
            }
            width1 = width;
            width = fabs(sty - stx);
        }

    }
}



#endif
