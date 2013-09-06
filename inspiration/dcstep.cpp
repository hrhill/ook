#include <iostream>
#include <cmath>
#include <algorithm>
#include <tuple>

#include "point.h"

inline
std::pair<double , double>
cubic_minimizer_(const point& p1, const point& p2)
{
	const double theta = 3.0 * (p1.f - p2.f) / (p2.st - p1.st) + p1.df + p2.df;
	const double gamma = copysign(sqrt(std::pow(theta, 2) - p1.df * p2.df), p2.st - p1.st);
	const double p = gamma - p1.df + theta;
	const double q = 2.0 * gamma - p1.df + p2.df;
	return std::make_pair(p / q, gamma);
}


inline
std::pair<double , double>
cubic_minimizer(double f1, double f2, double d1, double d2, double st1, double st2)
{
	/*
	const double theta = 3.0 * (f1 - f2) / (st2 - st1) + d1 + d2;
	const double gamma = copysign(sqrt(std::pow(theta, 2) - d1 * d2), st2 - st1);
	const double p = gamma - d1 + theta;
	const double q = 2.0 * gamma - d1 + d2;
	return std::make_pair(p / q, gamma);
	*/
	return cubic_minimizer_({f1, d1, st1}, {f2, d2, st2});
}

inline
std::pair<double, bool>
case1(double fx, double dx, double stx, double fp, double dp, double stp)
{
	double stpf, r, gamma;
	//std::tie(r, gamma) = cubic_minimizer(fx, fp, dx, dp, stx, stp);
	std::tie(r, gamma) = cubic_minimizer_({fx, dx, stx}, {fp, dp, stp});
	const double stpc = stx + r * (stp - stx);
	const double stpq = stx + dx / ((fx - fp) / (stp - stx) + dx) / 2. * (stp - stx);
	if (fabs(stpc - stx) < fabs(stpq - stx)){
	    stpf = stpc;
	} else {
	    stpf = stpc + (stpq - stpc) / 2.;
	}
	return std::make_pair(stpf, true);
}

inline
std::pair<double, bool>
case2(double fx, double dx, double stx, double fp, double dp, double stp)
{
	double stpf, r, gamma;
	//std::tie(r, gamma) = cubic_minimizer(fp, fx, dp, dx, stp, stx);
	std::tie(r, gamma) = cubic_minimizer_({fp, dp, stp}, {fx, dx, stx});	
	const double stpc = stp + r * (stx - stp);
	const double stpq = stp + dp / (dp - dx) * (stx - stp);
	if (fabs(stpc - stp) > fabs(stpq - stp)){
		stpf = stpc;
	} else {
		stpf = stpq;
	}
	return std::make_pair(stpf, true);	
}

inline
double
case3(double stx, double fx, double dx, double sty, double fy, double dy, double stp, double fp, double dp, bool brackt, double stpmin, double stpmax)
{
	/* 	The cubic step is computed only if the cubic tends to infinity
	in the direction of the step or if the minimum of the cubic
	is beyond stp. Otherwise the cubic step is defined to be the
	secant step. */
	double stpf, r, gamma;
	//std::tie(r, gamma) = cubic_minimizer(fp, fx, dp, dx, stp, stx);
	std::tie(r, gamma) = cubic_minimizer_({fp, dp, stp}, {fx, dx, stx});
	double stpc = (stp > stx) ? stpmax : stpmin;
	if (r < 0. && gamma != 0.) {
	    stpc = stp + r * (stx - stp);
	}
	const double delta = 0.66; // must be less than 1.	
	const double stpq = stp + dp / (dp - dx) * (stx - stp);
	if (brackt){
		/* A minimizer has been bracketed. If the cubic step is
		closer to stp than the secant step, the cubic step is
		taken, otherwise the secant step is taken. */
	    if (fabs(stpc - stp) < fabs(stpq - stp)) {
			stpf = stpc;
	    } else {
			stpf = stpq;
	    }
	    if (stp > stx) {
	    	stpf = std::min(stp + delta * (sty - stp), stpf);
	    } else {
			stpf = std::max(stp + delta * (sty - stp), stpf);
	    }
	} else {
		/* A minimizer has not been bracketed. If the cubic step is
		farther from stp than the secant step, the cubic step is
		taken, otherwise the secant step is taken. */
	    if (fabs(stpc - stp) > fabs(stpq - stp)) {
	    	stpf = stpc;
	    } else {
	    	stpf = stpq;
	    }
	    stpf = std::min(stpmax, stpf);
	    stpf = std::max(stpmin, stpf);
	}
	return stpf;
}

double
case4(double stx, double sty, double fy, double dy, double stp, double fp, double dp, bool brackt, double stpmin, double stpmax)
{
	double r, gamma;
	double stpf = (stp > stx) ? stpmax : stpmin;
	if (brackt) {
		//std::tie(r, gamma) = cubic_minimizer(fp, fy, dp, dy, stp, sty);
		std::tie(r, gamma) = cubic_minimizer_({fp, dp, stp}, {fy, dy, sty});		
	    stpf = stp + r * (sty - stp);
	}
	return stpf;
}

double
dcstep(double& stx, double& fx, double& dx, double& sty, double& fy, double& dy, const double& stp, const double& fp, const double& dp, bool& brackt, const double& stpmin, const double& stpmax)
{
/*	This subroutine computes a safeguarded step for a search
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

		subroutine dcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,stpmin,stpmax)

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
		On entry sty is the second endpoint of the interval that contains 
			the minimizer.
		On exit sty is the updated endpoint of the interval that contains 
			the minimizer.

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

	fp is the function at stp

	dp is the derivative of the function at stp.

	brackt is an bool variable.
		On entry brackt specifies if a minimizer has been bracketed.
			Initially brackt must be set to .false.
		On exit brackt specifies if a minimizer has been bracketed.
			When a minimizer is bracketed brackt is set to .true.

	stpmin is a lower bound for the step.

	stpmax is an upper bound for the step.
*/

    /* Local variables */
    double stpf;
    const double sgnd = dp * (dx / fabs(dx));
	/*	First case: A higher function value. The minimum is bracketed.
	If the cubic step is closer to stx than the quadratic step, the
	cubic step is taken, otherwise the average of the cubic and 
	quadratic steps is taken. */
    if (fp > fx) {
		std::tie(stpf, brackt) = case1(fx, dx, stx, fp, dp, stp);
	/*	Second case: A lower function value and derivatives of opposite
	sign. The minimum is bracketed. If the cubic step is farther from
	stp than the secant step, the cubic step is taken, otherwise the
	secant step is taken. */
    } else if (sgnd < 0.) {
		std::tie(stpf, brackt) = case2(fx, dx, stx, fp, dp, stp);
	/* Third case: A lower function value, derivatives of the same sign,
	and the magnitude of the derivative decreases. */
    } else if (fabs(dp) < fabs(dx)) {
		stpf = case3(stx, fx, dx, sty, fy, dy, stp, fp, dp, brackt, stpmin, stpmax);
	/* Fourth case: A lower function value, derivatives of the same sign,
	and the magnitude of the derivative does not decrease. If the
	minimum is not bracketed, the step is either stpmin or stpmax,
	otherwise the cubic step is taken. */
    } else {
    	stpf = case4(stx, sty, fy, dy, stp, fp, dp, brackt, stpmin, stpmax);
    }
	/* Update the interval which contains a minimizer. */
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
	/*     Compute the new step. */
    return stpf;
}