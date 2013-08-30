#include <cmath>
#include <algorithm>

int
dcstep_(double& stx, double& fx, double& dx, double& sty, double& fy, double& dy, double& stp, const double& fp, const double& dp, bool& brackt, double& stpmin, double& stpmax)
{
/*	Subroutine dcstep
	This subroutine computes a safeguarded step for a search
	procedure and updates an interval that contains a step that
	satisfies a sufficient decrease and a curvature condition.
	The parameter stx contains the step with the least function 
	value. If brackt is set to .true. then a minimizer has 
	been bracketed in an interval with endpoints stx and sty. 
	The parameter stp contains the current step. 
	The subroutine assumes that if brackt is set to .true. then 

		min(stx,sty) < stp < max(stx,sty), 

	and that the derivative at stx is negative in the direction 
	of the step. 

	The subroutine statement is 

		subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt, stpmin,stpmax)

	where

		stx is a double precision variable.
		On entry stx is the best step obtained so far and is an
			endpoint of the interval that contains the minimizer.
		On exit stx is the updated best step.

		fx is a double precision variable.
		On entry fx is the function at stx.
		On exit fx is the function at stx.

		dx is a double precision variable.
		On entry dx is the derivative of the function at 
			stx. The derivative must be negative in the direction of
			the step, that is, dx and stp - stx must have opposite signs. 
		On exit dx is the derivative of the function at stx.

		sty is a double precision variable. 
		On entry sty is the second endpoint of the interval that contains the 
			minimizer. 
		On exit sty is the updated endpoint of the interval that contains the 
			minimizer.

		fy is a double precision variable.
		On entry fy is the function at sty.
		On exit fy is the function at sty.

		dy is a double precision variable.
		On entry dy is the derivative of the function at sty. 
		On exit dy is the derivative of the function at the exit sty.

		stp is a double precision variable.
		On entry stp is the current step. If brackt is set to .true.
			then on input stp must be between stx and sty.
		On exit stp is a new trial step.

		fp is a double precision variable.
		On entry fp is the function at stp
		On exit fp is unchanged.

		dp is a double precision variable.
		On entry dp is the the derivative of the function at stp.
		On exit dp is unchanged.

		brackt is an bool variable.
		On entry brackt specifies if a minimizer has been bracketed. 
			Initially brackt must be set to .false.
		On exit brackt specifies if a minimizer has been bracketed.
			When a minimizer is bracketed brackt is set to .true.

		stpmin is a double precision variable.
		On entry stpmin is a lower bound for the step.
		On exit stpmin is unchanged.

		stpmax is a double precision variable.
		On entry stpmax is an upper bound for the step.
		On exit stpmax is unchanged. 

	MINPACK-1 Project. June 1983 
	Argonne National Laboratory. 
	Jorge J. More' and David J. Thuente.

	MINPACK-2 Project. November 1993.
	Argonne National Laboratory and University of Minnesota.
	Brett M. Averick and Jorge J. More'. */

    /* Local variables */
    static double stpc, stpf, stpq, p, q, r;

    double sgnd = dp * (dx / fabs(dx));
/*	First case: A higher function value. The minimum is bracketed. 
	If the cubic step is closer to stx than the quadratic step, the 
	cubic step is taken, otherwise the average of the cubic and 
	quadratic steps is taken. */
    if (fp > fx) {
    	const double theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;	
		const double s = std::max({fabs(theta), fabs(dx), fabs(dp)});	
	    const double gamma = copysign(s * sqrt(std::pow(theta / s, 2) - dy / s * (dp / s)), stp - stx);
		const double p = gamma - dx + theta;
		const double q = gamma - dx + gamma + dp;
		const double r = p / q;
		const double stpc = stx + r * (stp - stx);
		const double stpq = stx + dx / ((fx - fp) / (stp - stx) + dx) / 2.0 * (stp - stx);

		if (fabs(stpc - stx) < abs(stpq - stx)){
			stpf = stpc;
		} else {
			stpf = stpc + (stpq - stpc) / 2.0;
		}
		brackt = true;
	/*	Second case: A lower function value and derivatives of opposite
		sign. The minimum is bracketed. If the cubic step is farther from
		stp than the secant step, the cubic step is taken, otherwise the
		secant step is taken. */
	} else if (sgnd < 0.) {
    	const double theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;	
		const double s = std::max({fabs(theta), fabs(dx), fabs(dp)});	
	    const double gamma = copysign(s * sqrt(std::pow(theta / s, 2) - dy / s * (dp / s)), stx - stp);
		const double p = gamma - dp + theta;
		const double q = gamma - dp + gamma + dx;
		const double r = p / q;
		const double stpc = stp + r * (stx - stp);
		const double stpq = stp + dp / (dp - dx) * (stx - stp);

		if (fabs(stpc - stp) > fabs(stpq - stp)){
	    	stpf = stpc;
		} else {
	    	stpf = stpq;
		}
		brackt = true;
	/*	Third case: A lower function value, derivatives of the same sign,
		and the magnitude of the derivative decreases. */
    } else if (fabs(dp) < fabs(dx)) {
	/*	The cubic step is computed only if the cubic tends to infinity
		in the direction of the step or if the minimum of the cubic
		is beyond stp. Otherwise the cubic step is defined to be the
		secant step. 
		
		The case gamma = 0 only arises if the cubic does not tend 
		to infinity in the direction of the step. */
    	const double theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;	
		const double s = std::max({fabs(theta), fabs(dx), fabs(dp)});	
	    const double gamma = copysign(std::max(s * sqrt(std::pow(theta / s, 2) - dy / s * (dp / s)), 0.0), stx - stp);

		const double p = gamma - dp + theta;
		const double q = gamma + (dx - dp) + gamma;
		const double r = p / q;
		if (r < 0.0 && gamma != 0.) {
	    	stpc = stp + r * (stx - stp);
		} else if (stp > stx) {
	    	stpc = stpmax;
		} else {
	    	stpc = stpmin;
		}
		stpq = stp + dp / (dp - dx) * (stx - stp);
		if (brackt) {
		/*  A minimizer has been bracketed. If the cubic step is
			closer to stp than the secant step, the cubic step is
			taken, otherwise the secant step is taken. */
	    	if (fabs(stpc - stp) < fabs(stpq - stp)) {
				stpf = stpc;
	    	} else {
				stpf = stpq;
	    	}
	    	if (stp > stx) {
				stpf = std::min(stp + 2.0 / 3.0 * (sty - stp), stpf);
	    	} else {
				stpf = std::max(stp + 2.0 / 3.0 * (sty - stp), stpf);
	    	}
	    } else {
		/*	A minimizer has not been bracketed. If the cubic step is
			farther from stp than the secant step, the cubic step is
			taken, otherwise the secant step is taken. */
	    	if (fabs(stpc - stp) >  fabs(stpq - stp)){
				stpf = stpc;
	    	} else {
				stpf = stpq;
	    	}
	    	stpf = std::min(stpmax, stpf);
	    	stpf = std::max(stpmin, stpf);
		}
	/*	Fourth case: A lower function value, derivatives of the same sign,
		and the magnitude of the derivative does not decrease. If the
		minimum is not bracketed, the step is either stpmin or stpmax,
		otherwise the cubic step is taken. */
    } else {
		if (brackt) {
    		const double theta = 3.0 * (fp - fy) / (sty - stp) + dy + dp;
			const double s = std::max({fabs(theta), fabs(dx), fabs(dp)});	    		
	    	const double gamma = copysign(s * sqrt(std::pow(theta / s, 2) - dy / s * (dp / s)), sty - stp);
	    	const double p = gamma - dp + theta;
	    	const double q = gamma - dp + gamma + dy;
	    	const double r = p / q;
	    	stpc = stp + r * (sty - stp);
	    	stpf = stpc;
		} else if (stp > stx) {
	    	stpf = stpmax;
		} else {
	    	stpf = stpmin;
		}
    }
	/*     Update the interval which contains a minimizer. */
    if (fp > fx) {
		sty = stp;
		fy = fp;
		dy = dp;
    } else {
    	if (sgnd < 0.) {
    		sty = stx;
	    	fy = fx;
	    	dy = dx;
		}
		stx = stp;
		fx = fp;
		dx = dp;
    }
	/*     Compute the new step. */
    stp = stpf;
}

