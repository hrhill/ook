
// Copyright 2013 Harry Hill
//
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

#ifndef OOK_TEST_FUNCTIONS_PARABOLA_HPP_
#define OOK_TEST_FUNCTIONS_PARABOLA_HPP_

#include <limits>
#include <tuple>
#include <vector>

#include "ook/vector.hpp"

namespace ook
{
namespace test_functions
{

struct parabola
{

    std::tuple<double, vector, matrix>
    operator()(const vector& x) const
    {
        double f = 0.5 * std::pow(norm_2(x), 2);
        matrix id(n, n, 0);
        for (size_t i = 0; i < n; ++i)
            id(i, i) = 1.0;
        return std::make_tuple(f, x, id);
    }

    static const int n = 4;
    static const int m = 4;
    static double f_min;
    static double tolerance;
    static std::vector<double> local_minima;
    static std::vector<double> minima;
    static std::vector<double> x0;
};

double parabola::f_min = 0.0;
double parabola::tolerance = std::numeric_limits<double>::epsilon();

std::vector<double> parabola::minima = {0, 0, 0, 0};
std::vector<double> parabola::local_minima = {0, 0, 0, 0};
std::vector<double> parabola::x0 = {10.0, -4.1513, 1.0, 1000.0};

} // ns test_functions
} // ns ook

#endif
