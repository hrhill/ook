#ifndef OOK_LINE_SEARCH_SAFE_STEP_H_
#define OOK_LINE_SEARCH_SAFE_STEP_H_

#include <iostream>
#include <cmath>
#include <algorithm>
#include <tuple>

namespace ook{

namespace line_search {

template <typename T>
T
cubic_minimizer(T f1, T f2, T d1, T d2, T st1, T st2)
{
	const T theta = T(3.0) * (f1 - f2) / (st2 - st1) + d1 + d2;
	const T gamma = copysign(sqrt(std::pow(theta, 2) - d1 * d2), st2 - st1);
	const T p = gamma - d1 + theta;
	const T q = T(2.0) * gamma - d1 + d2;
	return p / q;
}

template <typename T>
inline
std::pair<T, bool>
case1(T fx, T dx, T stx, T fp, T dp, T stp)
{
	T stpf;
	const T r = cubic_minimizer(fx, fp, dx, dp, stx, stp);
	const T stpc = stx + r * (stp - stx);
	const T stpq = stx + dx / ((fx - fp) / (stp - stx) + dx) / T(2.0) * (stp - stx);
	if (fabs(stpc - stx) < fabs(stpq - stx)){
	    stpf = stpc;
	} else {
	    stpf = stpc + T(0.5) * (stpq - stpc);
	}
	return std::make_pair(stpf, true);
}

template <typename T>
inline
std::pair<T, bool>
case2(T fx, T dx, T stx, T fp, T dp, T stp)
{
	T stpf;
	const T r = cubic_minimizer(fp, fx, dp, dx, stp, stx);
	const T stpc = stp + r * (stx - stp);
	const T stpq = stp + dp / (dp - dx) * (stx - stp);
	if (fabs(stpc - stp) > fabs(stpq - stp)){
		stpf = stpc;
	} else {
		stpf = stpq;
	}
	return std::make_pair(stpf, true);	
}

template <typename T>
inline
T
case3(T stx, T fx, T dx, T sty, T fy, T dy, T stp, T fp, T dp, bool brackt, T stpmin, T stpmax)
{
	/* 	The cubic step is computed only if the cubic tends to infinity
	in the direction of the step or if the minimum of the cubic
	is beyond stp. Otherwise the cubic step is defined to be the
	secant step. */
	const T theta = 3.0 * (fx - fp) / (stp - stx) + dx + dp;
	/* The case gamma = 0 only arises if the cubic does not tend
	to infinity in the direction of the step.*/
	const T gamma = copysign(std::max(T(0.0), sqrt(std::pow(theta, 2) - dx * dp)), stx - stp);
	const T p = gamma - dp + theta;
	const T q = T(2.0) * gamma - dp + dx;
	const T r = p / q;
	
	T stpc = (stp > stx) ? stpmax : stpmin;
	if (r < 0. && gamma != T(0.0)) {
	    stpc = stp + r * (stx - stp);
	}
	T stpf;
	const T delta(0.66); // must be less than 1.	
	const T stpq = stp + dp / (dp - dx) * (stx - stp);
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

template <typename T>
T
case4(T stx, T sty, T fy, T dy, T stp, T fp, T dp, bool brackt, T stpmin, T stpmax)
{
	T stpf = (stp > stx) ? stpmax : stpmin;
	if (brackt) {
	    const T r = cubic_minimizer(fp, fy, dp, dy, stp, sty);
	    stpf = stp + r * (sty - stp);
	}
	return stpf;
}

template <typename T>
T
safe_step(T& stx, T& fx, T& dx, T& sty, T& fy, T& dy, const T& stp, const T& fp, const T& dp, bool& brackt, const T& stpmin, const T& stpmax)
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

	stx is a T precision variable.
		On entry stx is the best step obtained so far and is an
			endpoint of the interval that contains the minimizer.
		On exit stx is the updated best step.

	fx is a T precision variable.
		On entry fx is the function at stx.
		On exit fx is the function at stx.

	dx is a T precision variable.
		On entry dx is the derivative of the function at
			stx. The derivative must be negative in the direction of
			the step, that is, dx and stp - stx must have opposite signs.
		On exit dx is the derivative of the function at stx.

	sty is a T precision variable.
		On entry sty is the second endpoint of the interval that contains 
			the minimizer.
		On exit sty is the updated endpoint of the interval that contains 
			the minimizer.

	fy is a T precision variable.
		On entry fy is the function at sty.
		On exit fy is the function at sty.

	dy is a T precision variable.
		On entry dy is the derivative of the function at sty.
		On exit dy is the derivative of the function at the exit sty.

	stp is a T precision variable.
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
    T stpf;
    const T sgnd = dp * (dx / fabs(dx));
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
		if (sgnd < T(0.0)) {
	    	sty = stx;
	    	fy = fx;
	    	dy = dx;
		}
		stx = stp;
		fx = fp;
		dx = dp;
    }
    return stpf;
}

} //ns line_search

} // ns ook

#endif