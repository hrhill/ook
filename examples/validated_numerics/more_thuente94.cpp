#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iterator>

#include "ook/line_search/more_thuente.h"
#include "ook/options.h"
#include "ook/test_functions/line_search.h"
#include "ook/message.h"

using namespace std;

struct problem{
    int nprob;
    int ntries;
    double stp0;
    double ftol;
    double gtol;
    double xtol;
};

struct algo_result{
    ook::message info;
    int nfev;
    int ntry;
    int nprob;
    double stp0;
    double stp;
    double f;
    double g;

    friend
    std::ostream&
    operator<<(std::ostream& out, const algo_result& r)
    {
        out << setw(6) << r.nprob
            << setw(6) << r.ntry
            << setw(6) << r.nfev
            << setw(6) << (r.info == ook::message::convergence)
            << setw(11) << scientific << r.stp0
            << setw(10) << scientific << r.stp
            << setw(10) << scientific << r.f
            << setw(10) << scientific << r.g;
        return out;
    }
};

problem probs[] = {
    {1,    4,     1.e-3,    1.e-03,    1.e-01,    1.e-10},
    {2,    4,     1.e-3,    1.e-01,    1.e-01,    1.e-10},
    {3,    4,     1.e-3,    1.e-01,    1.e-01,    1.e-10},
    {4,    4,     1.e-3,    1.e-03,    1.e-03,    1.e-10},
    {5,    4,     1.e-3,    1.e-03,    1.e-03,    1.e-10},
    {6,    4,     1.e-3,    1.e-03,    1.e-03,    1.e-10}};

int main()
{
    cout.precision(2);
    for (const auto& p : probs){
        // Read in problem parameters.
        const int nprob = p.nprob;
        const int ntries = p.ntries;
        const double stp0 = p.stp0;
        const double ftol = p.ftol;
        const double gtol = p.gtol;
        const double xtol = p.xtol;

        double factor = 1.0;
        double g0, f0;
        vector<algo_result> results;

        for (int ntry = 0; ntry < ntries; ++ntry) {
            // Initialize the search.
            std::tie(f0, g0) = ook::test_functions::mtfcn(0.0, nprob);
            int nfev = 0;
            auto phi_trace = [&nfev, nprob](const double& x){
                ++nfev;
                return ook::test_functions::mtfcn(x, nprob);
            };

            double stp = factor * stp0;
            ook::options<double> opts(ftol, gtol, xtol, 0, 4.0 * std::max(1.0, stp));
            ook::message msg;
            double f, g;
            std::tie(msg, stp, f, g) = ook::line_search::more_thuente(phi_trace, f0, g0, stp, opts);

            // Record information on the algorithm.
            results.push_back({msg, nfev, ntry + 1, nprob, factor * stp0, stp, f, g});
            factor *= 100.;
        }
        cout << "\n\n Summary of  " << ntries << " calls to dcsrch\n"
                "\n  xtol = " <<  scientific << xtol <<
                "   ftol = " <<  scientific << ftol <<
                "   gtol = " <<  scientific << gtol <<
                "   g0 = " <<  scientific << g0 << "\n\n";
        cout << setw(7)  << "nprob" << setw(6) << "ntry"
             << setw(6)  << "nfev"  << setw(6) << "info"
             << setw(7)  << "x0"    << setw(10) << "x"
             << setw(10) << "f"     << setw(12) << "g\n\n";

        for (const auto& x : results) std::cout << x << "\n";
    }
    return 0;
}
