#ifndef OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_FREUDENSTEIN_ROTH_HPP_
#define OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_FREUDENSTEIN_ROTH_HPP_

#include <tuple>
#include <limits>
#include <vector>

#include "ook/vector.hpp"
#include "ook/matrix.hpp"

namespace ook{
namespace test_functions{

struct freudenstein_roth
{
    std::tuple<double, vector, matrix>
    operator()(const vector& x) const
    {
        vector df(2, 1.0);
        matrix d2f(2, 2, 0.0);

        const double x1 = x[0];
        const double x2 = x[1];

        const double f1 = -13 + x1 + ((5 - x2) * x2 - 2) * x2;
        const double f2 = -29 + x1 + ((x2 + 1) * x2 - 14) * x2;
        const double f = f1 * f1 + f2 * f2;

        const double df1x1 = 1;
        const double df1x2 = -3 * x2 * x2 + 10 * x2 - 2;
        const double d2f1x1 = 0;
        const double d2f1x2 = -6 * x2 + 10;
        const double d2f1x1x2 = 0;

        const double df2x1 = 1;
        const double df2x2 = 3 * x2 * x2 + 2 * x2 - 14;
        const double d2f2x1 = 0;
        const double d2f2x2 = 6 * x2 + 2;
        const double d2f2x1x2 = 0;

        df[0] = 2 * f1 * df1x1 + 2 * f2 * df2x1;
        df[1] = 2 * f1 * df1x2 + 2 * f2 * df2x2;

        d2f(0, 0) = 2 * f1 * d2f1x1 + 2 * df1x1 * df1x1 + 2 * f2 * d2f2x1 + 2 * df2x1 * df2x1;
        d2f(1, 1) = 2 * f1 * d2f1x2 + 2 * df1x2 * df1x2 + 2 * f2 * d2f2x2 + 2 * df2x2 * df2x2;

        d2f(0, 1) = d2f(1, 0) = 2 * df1x1 * df1x2 + 2 * f1 * d2f1x1x2 + 2 * df2x1 * df2x2 + 2 * f2 * d2f2x1x2;

        return std::make_tuple(f, df, d2f);
    }

    static const int n = 2;
    static const int m = 2;
    static double f_min;
    static double tolerance;
    static std::vector<double> minima;
    static std::vector<double> x0;
    static std::vector<double> local_minima;
};

double freudenstein_roth::f_min = 0.0;
double freudenstein_roth::tolerance = std::numeric_limits<double>::epsilon();
std::vector<double> freudenstein_roth::minima = {5.0, 4.0};
std::vector<double> freudenstein_roth::local_minima = {11.41278, -0.8968053};
std::vector<double> freudenstein_roth::x0 = {0.5,  2.0};

} // ns test_functions
} // ns ook

#endif
