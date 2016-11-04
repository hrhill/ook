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

#include <vector>
#include <tuple>

#include "ook/lbfgs.hpp"

std::tuple<double, std::vector<double>>
obj_fun(const std::vector<double>& x)
{
    const int n = x.size();
    std::vector<double> g(n);
    double f = 0.;
    for (int j = 0; j < n; j += 2){
        double t1 = 1.0 - x[j];
        double t2 = (x[j + 1] - x[j] * x[j]) * 10.;
        g[j + 1] = t2 * 20.0;
        g[j] = (x[j] * g[j + 1] + t1) * -2.0;
        f += t1 * t1 + t2 * t2;
    }
    return std::make_tuple(f, g);
}

int main()
{
    ook::lbfgs_options<double>
    opts(5, 1e-05, false, {{1, 0}}, 1e-04, 0.9, 2000, 20);

    int n = 100;
    std::vector<double> x(n);
    for (int j = 0; j < n; j += 2){
        x[j] = -1.2;
        x[j + 1] = 1.0;
    }

    auto soln = ook::lbfgs(obj_fun, x, opts);
}
