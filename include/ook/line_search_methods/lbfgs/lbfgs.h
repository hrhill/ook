#ifndef LBFGS_H_
#define LBFGS_H_

#include <algorithm>
#include <stdexcept>

#include <cstdio>
#include <cmath>
#include <cassert>

#include "mcstep.h"
#include "mcsrch.h"

extern "C"{
#include <cblas.h>
}

namespace ook{
namespace detail{

int
check_diag(double* diag, int n)
{
    for (int i = 0; i < n; ++i) {
        if (diag[i] <= 0.0) {
            throw std::runtime_error(" IFLAG= -2\n THE " + std::to_string(i) + "-TH DIAGONAL ELEMENT OF THE\n INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE");
        }
    }
    return 0;    
}

// label output, crap name change it.
int label(int* iprint, const int iter, const int nfun, const double gnorm, const int n, const int m, double* x, double f, double* g, const double stp, const bool finish)
{
    char fmt_20[] = "  N=%*d   NUMBER OF CORRECTIONS=%*d\n       INITIAL VALUES\n";
    char fmt_30[] = " F= %*.3E   GNORM= %*.3E\n";
    char fmt_70[] = "\n   I   NFN    FUNC        GNORM       STEPLENGTH\n";

    if (iter == 0) {
        puts("*************************************************");
        printf(fmt_20, 5, n, 2, m);
        printf(fmt_30, 10, f, 10, gnorm);
        puts("*************************************************");
        puts(fmt_70);
    } else {
        if (iprint[0] == 0 && (iter != 1 && ! finish)) {
            return 0;
        }
        if (iprint[0] != 0) {
            if ((iter - 1) % iprint[0] == 0 || finish) {
                if (iprint[1] > 1 && iter > 1) {
                    puts(fmt_70);
                }
                printf("%*d %*d %*.3E %*.3E %*.3E\n", 4, iter, 4, nfun, 13, f, 11, gnorm, 11, stp);                
            } else {
                return 0;
            }
        } else {
            if (iprint[1] > 1 && finish) {
                puts(fmt_70);
            }
            printf("%*d %*d %*.3E %*.3E %*.3E\n", 4, iter, 4, nfun, 13, f, 11, gnorm, 11, stp);                            
        }
        if (finish) {
            puts("\n THE MINIMIZATION TERMINATED WITHOUT DETECTING ERRORS.\n IFLAG = 0");
        }
    }
    return 0;
}

void
compute_hg(int& cp, int& point, double ys, int n, int m, double* w, double* g, int bound, double* diag, int ispt, int iypt)
{
    /*     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980, */
    /*     "Updating quasi-Newton matrices with limited storage", */
    /*     Mathematics of Computation, Vol.24, No.151, pp. 773-782. */
    cp = (point == 0) ? m : point;
    w[n + cp - 1] = 1.0 / ys;
    for (int i = 0; i < n; ++i) {
        w[i] = -g[i];
    }
    cp = point;
    for (int i = 0; i < bound; ++i) {
        cp = (cp == 0) ? m -1 : cp - 1;
        const double sq = cblas_ddot(n, &w[ispt + cp * n], 1, w, 1);
        int inmc = n + m + cp;
        const int iycn = iypt + cp * n;
        w[inmc] = w[n + cp] * sq;
        cblas_daxpy(n, -w[inmc], &w[iycn], 1, w, 1);
    }
    for (int i = 0; i < n; ++i) {
        w[i] *= diag[i];
    }
    for (int i = 0; i < bound; ++i) {
        const double yr = cblas_ddot(n, &w[iypt + cp * n], 1, w, 1);
        double beta = w[n + cp] * yr;
        int inmc = n + m + cp;
        beta = w[inmc] - beta;
        const int iscn = ispt + cp * n;
        cblas_daxpy(n, beta, &w[iscn], 1, w, 1);
        cp = (cp == m - 1) ? 0 : cp + 1;
    }
    /*     STORE THE NEW SEARCH DIRECTION */
    std::copy(w, w + n, w + ispt + point * n);    
}

template <typename T>
int lbfgs(int n, int m, T *x, T f, T *g, bool diagco, T *diag, int *iprint, T eps, T xtol, T *w, int& iflag)
{
/*        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION */
/*                          JORGE NOCEDAL */
/*                        *** July 1990 *** */
/*     This subroutine solves the unconstrained minimization problem */

/*                      min F(x),    x= (x1,x2,...,xN), */

/*      using the limited memory BFGS method. The routine is especially */
/*      effective on problems involving a large number of variables. In */
/*      a typical iteration of this method an approximation Hk to the */
/*      inverse of the Hessian is obtained by applying M BFGS updates to */
/*      a diagonal matrix Hk0, using information from the previous M steps. */
/*      The user specifies the number M, which determines the amount of */
/*      storage required by the routine. The user may also provide the */
/*      diagonal matrices Hk0 if not satisfied with the default choice. */
/*      The algorithm is described in "On the limited memory BFGS method */
/*      for large scale optimization", by D. Liu and J. Nocedal, */
/*      Mathematical Programming B 45 (1989) 503-528. */

/*      The user is required to calculate the function value F and its */
/*      gradient G. In order to allow the user complete control over */
/*      these computations, reverse  communication is used. The routine */
/*      must be called repeatedly under the control of the parameter */
/*      IFLAG. */

/*      The steplength is determined at each iteration by means of the */
/*      line search routine MCVSRCH, which is a slight modification of */
/*      the routine CSRCH written by More' and Thuente. */

/*      The calling statement is */

/*          CALL LBFGS(N,M,X,F,G,DIAGCO,DIAG,IPRINT,EPS,XTOL,W,IFLAG) */

/*      where */

/*     N       is an int variable that must be set by the user to the */
/*             number of variables. It is not altered by the routine. */
/*             Restriction: N>0. */

/*     M       is an int variable that must be set by the user to */
/*             the number of corrections used in the BFGS update. It */
/*             is not altered by the routine. Values of M less than 3 are */
/*             not recommended; large values of M will result in excessive */
/*             computing time. 3<= M <=7 is recommended. Restriction: M>0. */

/*     X       is a DOUBLE PRECISION array of length N. On initial entry */
/*             it must be set by the user to the values of the initial */
/*             estimate of the solution vector. On exit with IFLAG=0, it */
/*             contains the values of the variables at the best point */
/*             found (usually a solution). */

/*     F       is a DOUBLE PRECISION variable. Before initial entry and on */
/*             a re-entry with IFLAG=1, it must be set by the user to */
/*             contain the value of the function F at the point X. */

/*     G       is a DOUBLE PRECISION array of length N. Before initial */
/*             entry and on a re-entry with IFLAG=1, it must be set by */
/*             the user to contain the components of the gradient G at */
/*             the point X. */

/*     DIAGCO  is a bool variable that must be set to .TRUE. if the */
/*             user  wishes to provide the diagonal matrix Hk0 at each */
/*             iteration. Otherwise it should be set to .FALSE., in which */
/*             case  LBFGS will use a default value described below. If */
/*             DIAGCO is set to .TRUE. the routine will return at each */
/*             iteration of the algorithm with IFLAG=2, and the diagonal */
/*              matrix Hk0  must be provided in the array DIAG. */


/*     DIAG    is a DOUBLE PRECISION array of length N. If DIAGCO=.TRUE., */
/*             then on initial entry or on re-entry with IFLAG=2, DIAG */
/*             it must be set by the user to contain the values of the */
/*             diagonal matrix Hk0.  Restriction: all elements of DIAG */
/*             must be positive. */

/*     IPRINT  is an int array of length two which must be set by the */
/*             user. */

/*             IPRINT(0) specifies the frequency of the output: */
/*                IPRINT(0) < 0 : no output is generated, */
/*                IPRINT(0) = 0 : output only at first and last iteration, */
/*                IPRINT(0) > 0 : output every IPRINT(1) iterations. */

/*             IPRINT(2) specifies the type of output generated: */
/*                IPRINT(2) = 0 : iteration count, number of function */
/*                                evaluations, function value, norm of the */
/*                                gradient, and steplength, */
/*                IPRINT(2) = 1 : same as IPRINT(2)=0, plus vector of */
/*                                variables and  gradient vector at the */
/*                                initial point, */

/*     EPS     is a positive DOUBLE PRECISION variable that must be set by */
/*             the user, and determines the accuracy with which the solution */
/*             is to be found. The subroutine terminates when */

/*                         ||G|| < EPS std::max(1,||X||), */

/*             where ||.|| denotes the Euclidean norm. */

/*     XTOL    is a  positive DOUBLE PRECISION variable that must be set by */
/*             the user to an estimate of the machine precision (e.g. */
/*             10**(-16) on a SUN station 3/60). The line search routine will */
/*             terminate if the relative width of the interval of uncertainty */
/*             is less than XTOL. */

/*     W       is a DOUBLE PRECISION array of length N(2M+1)+2M used as */
/*             workspace for LBFGS. This array must not be altered by the */
/*             user. */

/*     IFLAG   is an int variable that must be set to 0 on initial entry */
/*             to the subroutine. A return with IFLAG<0 indicates an error, */
/*             and IFLAG=0 indicates that the routine has terminated without */
/*             detecting errors. On a return with IFLAG=1, the user must */
/*             evaluate the function F and gradient G. On a return with */
/*             IFLAG=2, the user must provide the diagonal matrix Hk0. */

/*             The following negative values of IFLAG, detecting an error, */
/*             are possible: */

/*              IFLAG=-1  The line search routine MCSRCH failed. The */
/*                        parameter INFO provides more detailed information */
/*                        (see also the documentation of MCSRCH): */

/*                       INFO = 0  IMPROPER INPUT PARAMETERS. */

/*                       INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF */
/*                                 UNCERTAINTY IS AT MOST XTOL. */

/*                       INFO = 3  MORE THAN 20 FUNCTION EVALUATIONS WERE */
/*                                 REQUIRED AT THE PRESENT ITERATION. */

/*                       INFO = 4  THE STEP IS TOO SMALL. */

/*                       INFO = 5  THE STEP IS TOO LARGE. */

/*                       INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS. */
/*                                 THERE MAY NOT BE A STEP WHICH SATISFIES */
/*                                 THE SUFFICIENT DECREASE AND CURVATURE */
/*                                 CONDITIONS. TOLERANCES MAY BE TOO SMALL. */

/*    ON THE DRIVER: */

/*     The subroutine contains one common area, which the user may wish to */
/*    reference: */

/*    MP  is an int variable with default value 6. It is used as the */
/*        unit number for the printing of the monitoring information */
/*        controlled by IPRINT. */

/*    LP  is an int variable with default value 6. It is used as the */
/*        unit number for the printing of error messages. This printing */
/*        may be suppressed by setting LP to a non-positive value. */

/*    GTOL is a DOUBLE PRECISION variable with default value 0.9, which */
/*        controls the accuracy of the line search routine MCSRCH. If the */
/*        function and gradient evaluations are inexpensive with respect */
/*        to the cost of the iteration (which is sometimes the case when */
/*        solving very large problems) it may be advantageous to set GTOL */
/*        to a small value. A typical small value is 0.1.  Restriction: */
/*        GTOL should be greater than 1.D-04. */

/*    STPMIN and STPMAX are non-negative DOUBLE PRECISION variables which */
/*        specify lower and uper bounds for the step in the line search. */
/*        Their default values are 1.D-20 and 1.D+20, respectively. These */
/*        values need not be modified unless the exponents are too large */
/*        for the machine being used, or unless the problem is extremely */
/*        badly scaled (in which case the exponents should be increased). */

/*    Other routines called directly:  MCSRCH */
    assert(n > 0 && m > 0);
    /* Local variables */
    static int info, nfev;
    static int ispt, iypt, bound;
    static int point;
    static int cp;
    static double ys;
    static bool finish;
    static double stp, stp1;

    /* Genuine static variables */
    static int iter = 0;
    static int nfun = 0;
    /* Local variables
       Parameters for line search, should be const */
    double ftol = 1e-4;
    int maxfev = 20;
    double gnorm = 0;
    int npt = 0;

    if (iflag == 0){
        nfun = 1;
        point = 0;
        finish = false;
        if (diagco) {
            iflag = check_diag(diag, n);
        } else {
            std::fill(diag, diag + n, 1.0);
        }
        /*     THE WORK VECTOR W IS DIVIDED AS FOLLOWS: */
        /*     --------------------------------------- */
        /*     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND */
        /*         OTHER TEMPORARY INFORMATION. */
        /*     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO. */
        /*     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED */
        /*         IN THE FORMULA THAT COMPUTES H*G. */
        /*     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH */
        /*         STEPS. */
        /*     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M */
        /*         GRADIENT DIFFERENCES. */

        /*     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A */
        /*     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT. */
        ispt = n + (m << 1);
        iypt = ispt + n * m;
        for (int i = 0; i < n; ++i) {
            w[ispt + i] = -g[i] * diag[i];
        }
        gnorm = sqrt(cblas_ddot(n, g, 1, g, 1));
        stp1 = 1.0 / gnorm;

        if (iprint[0] >= 0) {
            label(iprint, iter, nfun, gnorm, n, m, x, f, g, stp, finish);            
        }
    }else if (iflag == 1){
        mcsrch(n, x, f, g, &w[ispt + point * n], stp, ftol, xtol, maxfev, info, nfev, diag);
        if (info == -1) {
            return 0;
        }
        nfun += nfev;
        /*     COMPUTE THE NEW STEP AND GRADIENT CHANGE */
        npt = point * n;
        for (int i = 0; i < n; ++i) {
            w[ispt + npt + i] = stp * w[ispt + npt + i];
            w[iypt + npt + i] = g[i] - w[i];
        }
        point = (point == m - 1) ? 0 : point + 1;

        /* Termination test */
        gnorm = sqrt(cblas_ddot(n, g, 1, g, 1));
        double xnorm = sqrt(cblas_ddot(n, x, 1, x, 1));
        xnorm = std::max(1.0, xnorm);
        finish = (gnorm / xnorm <= eps);

        if (iprint[0] >= 0) {
            label(iprint, iter, nfun, gnorm, n, m, x, f, g, stp, finish);
        }
        if (finish) {
            iflag = 0;
            return 0;
        }
    }

    while(1){
        info = 0;        
        if (iflag == 0 || iflag == 1){
            ++iter;
            bound = (iter > m) ? m : iter - 1;        
            if (iter != 1) {
                ys = cblas_ddot(n, &w[iypt + npt], 1, &w[ispt + npt], 1);
                if (diagco){
                    iflag = 2;
                    return 0;                
                }else{
                    const double yy = cblas_ddot(n, &w[iypt + npt], 1, &w[iypt + npt], 1);
                    std::fill(diag, diag + n, ys/yy);
                }
                compute_hg(cp, point, ys, n, m, w, g, bound, diag, ispt, iypt);
            }
        }else if (iflag == 2){
            iflag = check_diag(diag, n);            
            compute_hg(cp, point, ys, n, m, w, g, bound, diag, ispt, iypt);
        }
        /*     OBTAIN THE ONE-DIMENSIONAL MINIMIZER OF THE FUNCTION */
        /*     BY USING THE LINE SEARCH ROUTINE MCSRCH */
        nfev = 0;
        stp = (iter == 1) ? stp1 : 1.0;
        std::copy(g, g + n, w);
        mcsrch(n, x, f, g, &w[ispt + point * n], stp, ftol, xtol, maxfev, info, nfev, diag);
        if (info == -1) {
            iflag = 1;
            return 0;
        }
        nfun += nfev;

        /*     COMPUTE THE NEW STEP AND GRADIENT CHANGE */
        npt = point * n;
        for (int i = 0; i < n; ++i) {
            w[ispt + npt + i] = stp * w[ispt + npt + i];
            w[iypt + npt + i] = g[i] - w[i];
        }
        point = (point == m - 1) ? 0 : point + 1;

        /* Termination test */
        gnorm = sqrt(cblas_ddot(n, g, 1, g, 1));
        double xnorm = sqrt(cblas_ddot(n, x, 1, x, 1));
        xnorm = std::max(1.0, xnorm);
        finish = (gnorm / xnorm <= eps);

        if (iprint[0] >= 0) {
            label(iprint, iter, nfun, gnorm, n, m, x, f, g, stp, finish);
        }
        if (finish) {
            iflag = 0;
            return 0;
        }
    }
}

} // ns detail
} // ns ook


#endif