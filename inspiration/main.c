#include <stdio.h>
#include <stdbool.h>

#include "line_search.h"

#define max(x, y) ((x) > (y) ? (x) : (y))
#define min(x, y) ((x) < (y) ? (x) : (y))

double
phi51(double a, double b, double* df){
    *df = (a * a - b)/((a * a + b) * (a * a + b));
    return -a/(a * a + b);
}

/*
double
phi52(double a, double b, double* df){
    const double apb = (a + b);

    return apb * (1.0 - 2.0 * apb * apb * apb * apb);
}
*/
int main(int argc, char** argv){

    const double b = 2.0;
    const double factor = 1.0;
    const int ntries = 100;

    double stp = 0.0;
    double g = 0;

    double f = phi51(stp, b, &g);
    double stpmin = 0.0;
    double stpmax = 4.0 * max(1.0, stp);
    double ftol = 1e-03;
    double gtol = 1e-01;
    double xtol = 1e-10;
    double stp0 = 1e-03;

    stp = factor * stp0;
    int nfev = 0;
    char task[255];
    s_copy(task, "START", 5, 5);
    int isave[2];
    double dsave[13];
    int i;
    for (i = 0; i < ntries; ++i){
        dcsrch_(&stp, &f, &g, &ftol, &gtol, &xtol, task, &stpmin, &stpmax, isave, dsave, 255);
        if (s_cmp(task, "FG", 2, 2) == 0){
            nfev = nfev + 1;
            f = phi51(stp, b, &g);
        }
        printf("%d, %f, %f, %f, %s\n", nfev, stp, f, g, task);
    };

    return 0;
}