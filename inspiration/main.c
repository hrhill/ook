#include <iostream>
#include <string>
#include <algorithm>
#include <tuple>

#include "line_search.h"

std::tuple<double, double>
phi51(double a, double b){
    return std::make_tuple(-a/(a * a + b),
                          (a * a - b)/((a * a + b) * (a * a + b)));
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
    double f, g;

    std::tie(f, g) = phi51(stp, b);

    double ftol = 1e-03;
    double gtol = 1e-01;
    double xtol = 1e-10;
    double stp0 = 1e-03;

    double stpmin = 0.0;
    double stpmax = 4.0 * std::max(1.0, stp0);

    stp = factor * stp0;
    int nfev = 0;
    std::string task("START");

    for (int i = 0; i < ntries; ++i){
        dcsrch_(stp, f, g, ftol, gtol, xtol, task, stpmin, stpmax);
        if (task.find("FG") != std::string::npos) {
            nfev = nfev + 1;
            std::tie(f, g) = phi51(stp, b);
        }
        std::cout << nfev << ", " <<  stp  << ", " << f  << ", " << g  << ", " <<  task << std::endl;
        if (task.find("CONV") != std::string::npos)
            break;
    };

    return 0;
}