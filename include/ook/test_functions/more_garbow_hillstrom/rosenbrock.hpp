#ifndef OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_ROSENBROCK_HPP_
#define OOK_TEST_FUNCTIONS_MORE_GARBOW_HILLSTROM_ROSENBROCK_HPP_

#include <tuple>
#include <limits>
#include <vector>

#include "ook/vector.hpp"
#include "ook/matrix.hpp"

namespace ook{
namespace test_functions{

struct rosenbrock
{
    std::tuple<double, vector, matrix>
    operator()(const vector& x) const
    {
        const double x1 = x[0];
        const double x2 = x[1];
        const double f1 = 10 * (x2 - x1 * x1);
        const double f2 = 1- x1;
        const double f = f1 * f1 + f2 * f2;

        vector df(2, 1.0);
        df[0] = - 2 * f1 * 20 * x1 - 2 * f2;
        df[1] = 2 * f1 * 10;

        matrix d2f(2, 2, 0.0);
        d2f(0, 0) = 2 - 400 * x2 + 1200 * x1 * x1;
        d2f(1, 0) = d2f(0, 1) = -400 * x1;
        d2f(1, 1) = 200;

        return std::make_tuple(f, df, d2f);
    }

    static const int n = 2;
    static const int m = 2;
    static double f_min;
    static double tolerance;
    static std::vector<double> minima;
    static std::vector<double> local_minima;
    static std::vector<double> x0;
};

double rosenbrock::f_min = 0.0;
double rosenbrock::tolerance = std::numeric_limits<double>::epsilon();
std::vector<double> rosenbrock::minima = {1.0, 1.0};
std::vector<double> rosenbrock::local_minima = {1.0, 1.0};
std::vector<double> rosenbrock::x0 = {-1.2, 1.0};


} // ns test_functions
} // ns ook

#endif
