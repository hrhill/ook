// This file is part of ook.
//
// ook is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ook is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with ook.  If not, see <http://www.gnu.org/licenses/>.

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <cmath>
#include <algorithm>
#include <vector>
#include <iterator>

#include "ook/line_search/mcsrch.hpp"
#include "ook/test_functions/line_search.hpp"
#include "ook/line_search/options.hpp"
#include "ook/message.hpp"

using namespace std;

struct problem{
    int nprob;
    int ntries;
    double stp0;
    double ftol;
    double gtol;
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
    {1,    4,     1.e-3,    1.e-03,    1.e-01},
    {2,    4,     1.e-3,    1.e-01,    1.e-01},
    {3,    4,     1.e-3,    1.e-01,    1.e-01},
    {4,    4,     1.e-3,    1.e-03,    1.e-03},
    {5,    4,     1.e-3,    1.e-03,    1.e-03},
    {6,    4,     1.e-3,    1.e-03,    1.e-03}};

int main()
{
    cout.precision(2);
    cout <<
        "\nThis program validates that the translation of the More Thuente\n"
        "line search algorithm remains faithful to the implementation described\n"
        "in the paper\n"
        "author = {Jorge J. More and David J. Thuente and Preprint Mcs-p},\n"
        "title = {Line Search Algorithms With Guaranteed Sufficient Decrease},\n"
        "journal = {ACM Trans. Math. Software},\n"
        "year = {1992},\n"
        "volume = {20},\n"
        "pages = {286--307}\n"
        "Press enter to continue\n\n";
    cin.get();

    cout.precision(2);
    vector<string> table{"Table 5.1, ftol = 0.001, gtol = 0.1\n",
                         "Table 5.2, ftol = 0.1,   gtol = 0.1\n",
                         "Table 5.3, ftol = 0.1,   gtol = 0.1\n",
                         "Table 5.4, ftol = 0.001, gtol = 0.001\n",
                         "Table 5.5, ftol = 0.001, gtol = 0.001\n",
                         "Table 5.6, ftol = 0.001, gtol = 0.001\n"};
    string line(80, '-');

    for (const auto& p : probs){
        // Read in problem parameters.
        const int nprob = p.nprob;
        const int ntries = p.ntries;
        const double stp0 = p.stp0;
        const double ftol = p.ftol;
        const double gtol = p.gtol;

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
            ook::line_search::options<double> opts(ftol, gtol, 0, 4.0 * std::max(1.0, stp));
            ook::message msg;
            double f, g;
            tie(msg, stp, f, g) = ook::line_search::mcsrch(phi_trace, f0, g0, stp, opts);

            // Record information on the algorithm.
            results.push_back({msg, nfev, ntry + 1, nprob, factor * stp0, stp, f, g});
            factor *= 100.;
        }
        cout << line << "\n" << table[p.nprob - 1] << endl;
        cout << " Summary of  " << ntries << " calls to dcsrch\n"
                "   ftol = " <<  scientific << ftol <<
                "   gtol = " <<  scientific << gtol <<
                "   g0 = " <<  scientific << g0 << "\n\n";
        cout << setw(7)  << "nprob" << setw(6) << "ntry"
             << setw(6)  << "nfev"  << setw(6) << "info"
             << setw(7)  << "x0"    << setw(10) << "x"
             << setw(10) << "f"     << setw(12) << "g\n\n";

        for (const auto& x : results) cout << x << "\n";
    }
    return 0;
}
