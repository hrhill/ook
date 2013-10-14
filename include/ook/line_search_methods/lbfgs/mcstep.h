#ifndef MCSTEP_H_
#define MCSTEP_H_

#include <stdexcept>

template <typename T>
int
mcstep(T& stx, T& fx, T& dx, T& sty, T& fy, T& dy, T& stp, T& fp, T& dp, bool& brackt, const T& stpmin, const T& stpmax)
{
    /*     THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR */
    /*     A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR */
    /*     A MINIMIZER OF THE FUNCTION. */

    /*     THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION */
    /*     VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS */
    /*     ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE */
    /*     DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A */
    /*     MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY */
    /*     WITH ENDPOINTS STX AND STY. */

    /*     THE SUBROUTINE STATEMENT IS */

    /*       SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,STPMIN,STPMAX) */

    /*     WHERE */

    /*       STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP, */
    /*         THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED */
    /*         SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION */
    /*         OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE */
    /*         SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY. */

    /*       STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP, */
    /*         THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF */
    /*         THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE */
    /*         UPDATED APPROPRIATELY. */

    /*       STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP, */
    /*         THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP. */
    /*         IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE */
    /*         BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP. */

    /*       BRACKT IS A bool VARIABLE WHICH SPECIFIES IF A MINIMIZER */
    /*         HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED */
    /*         THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER */
    /*         IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE. */

    /*       STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER */
    /*         AND UPPER BOUNDS FOR THE STEP. */

    /*     CHECK THE INPUT PARAMETERS FOR ERRORS. */
    if ((brackt && (stp <= std::min(stx,sty) || stp >= std::max(stx,sty))) || dx *(stp - stx) >= 0.0 || stpmax < stpmin) {
        throw std::runtime_error("invalid parameters passed to mcstep");
    }

    double stpf;
    bool bound;

    /*     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN. */
    double sgnd = dp * (dx / fabs(dx));

    /*     FIRST CASE. A HIGHER FUNCTION VALUE. */
    /*     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER */
    /*     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN, */
    /*     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN. */
    if (fp > fx) {
        double theta = (fx - fp) * 3 / (stp - stx) + dx + dp;
        double s = std::max({fabs(theta), fabs(dy), fabs(dp)});
        double gamma = s * sqrt(theta / s * theta / s - dx / s * (dp / s));
        if (stp < stx) {
            gamma = -gamma;
        }
        double p = gamma - dx + theta;
        double q = gamma - dx + gamma + dp;
        double r = p / q;
        double stpc = stx + r * (stp - stx);
        double stpq = stx + dx / ((fx - fp) / (stp - stx) + dx) / 2 * (stp - stx);
        if (fabs(stpc - stx) < fabs(stpq - stx)){
            stpf = stpc;
        } else {
            stpf = stpc + (stpq - stpc) / 2;
        }
        bound = true;        
        brackt = true;
    /*     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF */
    /*     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC */
    /*     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP, */
    /*     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN. */
    } else if (sgnd < 0.f) {
        double theta = (fx - fp) * 3 / (stp - stx) + dx + dp;
        double s = std::max({fabs(theta), fabs(dy), fabs(dp)});        
        double gamma = s * sqrt(theta / s * theta / s - dx / s * (dp / s));
        if (stp > stx) {
            gamma = -gamma;
        }
        double p = gamma - dp + theta;
        double q = gamma - dp + gamma + dx;
        double r = p / q;
        double stpc = stp + r * (stx - stp);
        double stpq = stp + dp / (dp - dx) * (stx - stp);
        if (fabs(stpc - stp) > fabs(stpq - stp)){
            stpf = stpc;
        } else {
            stpf = stpq;
        }
        bound = false;        
        brackt = true;
    /*     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE */
    /*     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES. */
    /*     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY */
    /*     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC */
    /*     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE */
    /*     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO */
    /*     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP */
    /*     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN. */
    } else if (fabs(dp) < fabs(dx)) {
        double theta = (fx - fp) * 3 / (stp - stx) + dx + dp;
        double s = std::max({fabs(theta), fabs(dy), fabs(dp)});
        /*        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND */
        /*        TO INFINITY IN THE DIRECTION OF THE STEP. */
        double gamma = s * sqrt(std::max(0.0, theta / s * theta / s - dx / s * (dp / s)));
        if (stp > stx) {
            gamma = -gamma;
        }
        double p = gamma - dp + theta;
        double q = gamma + (dx - dp) + gamma;
        double r = p / q;
        double stpc = stpmin;
        if (r < 0.0 && gamma != 0.0) {
            stpc = stp + r * (stx - stp);
        } else if (stp > stx) {
            stpc = stpmax;
        }
        double stpq = stp + dp / (dp - dx) * (stx - stp);
        if (brackt) {
            if (fabs(stp - stpc) < fabs(stp - stpq)) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
        } else {
            if (fabs(stp - stpc) > fabs(stp - stpq)) {
                stpf = stpc;
            } else {
                stpf = stpq;
            }
        }
        bound = true;
    /*     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE */
    /*     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES */
    /*     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP */
    /*     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN. */
    } else {
        bound = false;
        if (brackt) {
            double theta = (fp - fy) * 3 / (sty - stp) + dy + dp;
            double s = std::max({fabs(theta), fabs(dy), fabs(dp)});
            double gamma = s * sqrt(theta / s * theta / s - dy / s * (dp / s));
            if (stp > sty) {
                gamma = -gamma;
            }
            double p = gamma - dp + theta;
            double q = gamma - dp + gamma + dy;
            double r = p / q;
            double stpc = stp + r * (sty - stp);
            stpf = stpc;
        } else if (stp > stx) {
            stpf = stpmax;
        } else {
            stpf = stpmin;
        }
    }
    /*     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT */
    /*     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE. */
    if (fp > fx) {
        sty = stp;
        fy = fp;
        dy = dp;
    } else {
        if (sgnd < 0.0) {
            sty = stx;
            fy = fx;
            dy = dx;
        }
        stx = stp;
        fx = fp;
        dx = dp;
    }
    /*     COMPUTE THE NEW STEP AND SAFEGUARD IT. */
    stpf = std::min(stpmax,stpf);
    stpf = std::max(stpmin,stpf);
    stp = stpf;
    if (brackt && bound) {
        if (sty > stx) {
            stp = std::min(stx + (sty - stx) * .66f, stp);
        } else {
            stp = std::max(stx + (sty - stx) * .66f, stp);
        }
    }
    return 0;
}

#endif