#ifndef OOK_LINE_SEARCH_MORE_THUENTE_H_
#define OOK_LINE_SEARCH_MORE_THUENTE_H_

#include <tuple>
#include <exception>
#include <cassert>

#include "ook/message.h"
#include "ook/line_search/conditions.h"
#include "ook/line_search/more_thuente/search.h"

namespace ook{
namespace line_search{

struct more_thuente{

    template <typename F, typename T, typename Options>
    std::tuple<message, T, T, T>
    operator()(F phi, T phi0, T dphi0, T a, const Options& opts)
    {
        if(a < opts.stpmin)
            return std::make_tuple(message::warning_stp_eq_stpmin, 0, phi0, dphi0);
        if(a > opts.stpmax)
            return std::make_tuple(message::warning_stp_eq_stpmax, 0, phi0, dphi0);
        if (dphi0 > T(0.0)) {
            return std::make_tuple(message::search_direction_is_not_a_descent_direction, 0, phi0, dphi0);
        }
        /*
        bool brackt;
        int stage;
        double gx;
        double gy;
        double fx;
        double fy;
        double stx;
        double sty;
        double stmin;
        double stmax;
        double width;
        double width1;
        */
        T phia, dphia;
        std::tie(phia, dphia) = phi(a);
        detail::state s = {false, 1, dphi0, dphia, phi0, phia, 0, a, 0, 5 * a,
                            (opts.stpmax - opts.stpmin),
                            (opts.stpmax - opts.stpmin) / 0.5};
        while(true){
            // Test for convergence.
            const bool sufficient_decrease = sufficient_decrease_condition(phia, phi0, opts.ftol, a, dphi0);
            const bool curvature = curvature_condition(dphia, opts.gtol, dphi0);

            if (sufficient_decrease && curvature)
                break;
            // Test for warnings.
            if (s.brackt && (a <= s.stmin || a >= s.stmax)) {
                return std::make_tuple(message::warning_rounding_error_prevents_progress, a, phia, dphia);
            }
            if (s.brackt && s.stmax - s.stmin <= opts.xtol * s.stmax) {
                return std::make_tuple(message::warning_xtol_satisfied, a, phia, dphia);
            }
            if (a == opts.stpmax && sufficient_decrease && curvature) {
                return std::make_tuple(message::warning_stp_eq_stpmax, a, phia, dphia);
            }
            if (a == opts.stpmin && (!sufficient_decrease || !curvature)) {
                return std::make_tuple(message::warning_stp_eq_stpmin, a, phia, dphia);
            }
            a = detail::search(a, phia, dphia, phi0, dphi0, opts, s);
            std::tie(phia, dphia) = phi(a);
        }
        return std::make_tuple(message::convergence, a, phia, dphia);
    }
};

} // ns linesearch
} // ns ook

#endif
