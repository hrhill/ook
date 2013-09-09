#ifndef OOK_MORE_THUENTE_H_
#define OOK_MORE_THUENTE_H_

#include <tuple>
#include <iomanip>

#include "dcsrch.h"

struct
more_thuente_state{
  double stp;
  task_value task;
  int nfev;
  double f;
  double g;
};

template <typename F>
more_thuente_state
more_thuente_line_search(F phi, double stp, const options& opts){

    task_value task = task_value::start;
    double f0, g0;
    std::tie(f0, g0) = phi(0.0);

    int nfev = 1;

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

        std::tie(task, stp) = dcsrch_(stp, f, g, opts);
        ++nfev;
        std::tie(f, g) = phi(stp);
    }while (task == task_value::update);

    return more_thuente_state{stp, task, nfev, f, g};
}

#endif