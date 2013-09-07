#ifndef OOK_MORE_THUENTE_H_
#define OOK_MORE_THUENTE_H_

#include <tuple>
#include <iomanip>

#include "dcsrch.h"

template <typename F>
void
more_thuente_line_search(F phi, const double stp0, const options& opts){

    double stp = 0.0;
    double f0, g0;
    std::tie(f0, g0) = phi(stp);

    stp = stp0;
    int nfev = 1;
    task_value task = task_value::start;

    double f, g;
    std::tie(f, g) = phi(stp);
    dcsrch_struct dcsrch_(f0, g0, stp, opts.stpmax - opts.stpmin);
    do{
        const bool sufficient_decrease = sufficient_decrease_condition(f, f0, g0, stp, opts.ftol);
        const bool curvature = curvature_condition(g, g0, opts.gtol);
        if (sufficient_decrease && curvature){
            task = task_value::convergence;
            break;
        }
        task = dcsrch_(stp, f, g, opts);
        if (task == task_value::fg) {
            nfev = nfev + 1;
            std::tie(f, g) = phi(stp);
        }
    }while (true);
    std::cout << std::scientific 
              << std::setw(16) << stp0 
              << std::setw(16) << task 
              << std::setw(4) << nfev 
              << std::setw(16) << stp  
              << std::setw(16) << g << std::endl;    
}

#endif