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

#ifndef OOK_TEST_FUNCTIONS_PARABOLA_H_
#define OOK_TEST_FUNCTIONS_PARABOLA_H_

#include <tuple>
#include <limits>
#include <vector>

#include "ook/norms.h"

namespace ook{
namespace test_functions{

template <typename Vector>
struct parabola
{
    typedef Vector vector_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type>
    operator()(const vector_type& x) const
    {
        real_type f = 0.5 * std::pow(norm_2(x), 2);
        return std::make_pair(f, x);
    }

    static const int n = 4;
    static const int m = 4;
    static real_type f_min;
    static real_type tolerance;
    static std::vector<real_type> minima;
    static std::vector<real_type> x0;
};

template <typename Vector>
typename Vector::value_type
parabola<Vector>::f_min = 0.0;

template <typename Vector>
typename Vector::value_type
parabola<Vector>::tolerance = std::numeric_limits<typename Vector::value_type>::epsilon();

template <typename Vector>
std::vector<typename Vector::value_type>
parabola<Vector>::minima = {0, 0, 0, 0};

template <typename Vector>
std::vector<typename Vector::value_type>
parabola<Vector>::x0 = {10.0, -4.1513, 1.0, 1000.0};


} // ns test_functions
} // ns ook

#endif
