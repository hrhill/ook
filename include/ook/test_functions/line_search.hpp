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

#ifndef OOK_TEST_FUNCTIONS_LINE_SEARCH_FUNCTIONS_HPP_
#define OOK_TEST_FUNCTIONS_LINE_SEARCH_FUNCTIONS_HPP_

#include <boost/math/constants/constants.hpp>
#include <cmath>
#include <tuple>

namespace ook
{

namespace test_functions
{

template <typename T>
std::tuple<T, T>
constant(const T& x)
{
    return std::make_tuple(-1.0, 0.0);
}

template <typename T>
std::tuple<T, T>
linear(const T& x, const T& a, const T& b)
{
    return std::make_tuple(a * x + b, a);
}

template <typename T>
std::tuple<T, T>
phi51(const T& a, const T& b)
{
    return std::make_tuple(-a / (a * a + b),
                           (a * a - b) / ((a * a + b) * (a * a + b)));
}

template <typename T>
std::tuple<T, T>
phi52(const T& a, const T& b)
{
    T t = a + b;
    T t2 = std::pow(t, 2);
    T t3 = std::pow(t, 3);
    T t4 = std::pow(t2, 2);
    T f = t * t4 - 2 * t4;
    T g = 5 * t4 - 8 * t3;
    return std::make_tuple(f, g);
}

template <typename T>
std::tuple<T, T>
phi53(const T& a, const T& b, const T& c)
{
    const T pi = boost::math::constants::pi<T>();

    T phi0;
    T dphi0;
    if (a <= 1 - b)
    {
        phi0 = 1 - a;
        dphi0 = -1;
    }
    else if (a >= 1 + b)
    {
        phi0 = a - 1;
        dphi0 = 1;
    }
    else
    {
        phi0 = 0.5 * std::pow(a - 1, 2) / b + 0.5 * b;
        dphi0 = (a - 1) / b;
    }
    return std::make_tuple(phi0 +
                               2 * (1 - b) / (c * pi) * sin(0.5 * c * pi * a),
                           dphi0 + (1 - b) * cos(0.5 * c * pi * a));
}

template <typename T>
std::tuple<T, T>
phi54(const T& a, const T& b1, const T& b2)
{
    using std::pow;
    T t1 = sqrt(pow(b1, 2) + 1) - b1;
    T t2 = sqrt(pow(b2, 2) + 1) - b2;
    T k1 = sqrt(pow(1 - a, 2) + b2 * b2);
    T k2 = sqrt(pow(a, 2) + pow(b1, 2));
    T f = t1 * k1 + t2 * k2;
    T g = -t1 * (1 - a) / k1 + t2 * a / k2;
    return std::make_tuple(f, g);
}

template <typename T>
std::tuple<T, T>
mtfcn(const T& x, int nprob)
{
    std::tuple<T, T> f;
    if (nprob == 1)
    {
        f = phi51(x, 2.0);
    }
    else if (nprob == 2)
    {
        f = phi52(x, 0.004);
    }
    else if (nprob == 3)
    {
        f = phi53(x, 0.01, 39.0);
    }
    else if (nprob == 4)
    {
        f = phi54(x, 1e-3, 1e-3);
    }
    else if (nprob == 5)
    {
        f = phi54(x, 1e-2, 1e-3);
    }
    else if (nprob == 6)
    {
        f = phi54(x, 1e-3, 1e-2);
    }
    return f;
}

} // ns test_functions
} // ns ook

#endif
