
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

#include <tuple>
#include <limits>
#include <vector>

#include "linalg/norms.hpp"

namespace ook{
namespace test_functions{

template <typename Vector, typename Matrix>
struct parabola
{
    typedef Vector vector_type;
    typedef Matrix matrix_type;
    typedef typename vector_type::value_type real_type;

    std::tuple<real_type, vector_type, matrix_type>
    operator()(const vector_type& x) const
    {
        real_type f = 0.5 * std::pow(norm_2(x), 2);
        matrix_type id(n, n, 0);
        for (size_t i = 0; i < n; ++i) id(i, i) = 1.0;
        return std::make_tuple(f, x, id);
    }

    static const int n = 4;
    static const int m = 4;
    static real_type f_min;
    static real_type tolerance;
    static std::vector<real_type> local_minima;
    static std::vector<real_type> minima;
    static std::vector<real_type> x0;
};

template <typename Vector, typename Matrix>
typename Vector::value_type
parabola<Vector, Matrix>::f_min = 0.0;

template <typename Vector, typename Matrix>
typename Vector::value_type
parabola<Vector, Matrix>::tolerance = std::numeric_limits<typename Vector::value_type>::epsilon();

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type>
parabola<Vector, Matrix>::minima = {0, 0, 0, 0};

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type>
parabola<Vector, Matrix>::local_minima = {0, 0, 0, 0};

template <typename Vector, typename Matrix>
std::vector<typename Vector::value_type>
parabola<Vector, Matrix>::x0 = {10.0, -4.1513, 1.0, 1000.0};


} // ns test_functions
} // ns ook

#endif
