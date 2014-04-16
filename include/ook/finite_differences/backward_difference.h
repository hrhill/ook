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

#ifndef OOK_FINITE_DIFFERENCES_BACKWARD_DIFFERENCE_H_
#define OOK_FINITE_DIFFERENCES_BACKWARD_DIFFERENCE_H_

#include <limits>
#include <cmath>
#include <tuple>

#include <algorithm>
#include "ook/finite_differences/detail/transform.h"

namespace ook{
namespace finite_differences{

/// \brief Backward difference approximation.
struct backward_difference{
	/// \brief Calculate the backward difference approximation to the gradient of f.
	template <typename F, typename X>
	static
	std::tuple<typename X::value_type, X>
	gradient(F f, X x);

	/// \brief Calculate the backward difference approximation to the hessian of f.
	template <typename F, typename X, typename M>
	static
	std::tuple<typename X::value_type, M>
	hessian(F f, const X& x);
};

template <typename F, typename X>
std::tuple<typename X::value_type, X>
backward_difference::gradient(F f, X x)
{
	typedef typename X::value_type value_type;
	typedef typename X::size_type size_type;

	// Generate a set of sample points
	// (x + he_1, x + he_2, ..., x)
	const size_type n = std::distance(x.begin(), x.end());
	const value_type hmin(sqrt(std::numeric_limits<value_type>::epsilon()));
	std::vector<X> sample_points(n + 1);
	X h(n);

	for (size_type i = 0; i < n; ++i)
	{
		const value_type xi = x[i];
		const value_type hx = hmin * (1 + fabs(xi));
		h[i] = hx;
		x[i] = xi - hx;
		sample_points[i] = x;
		x[i] = xi;
	}
	sample_points[n] = x;

	// evaluate function at each point
	std::vector<value_type> function_values(sample_points.size());
	ook::finite_differences::detail::transform(sample_points.begin(), sample_points.end(), function_values.begin(), f);

	// assemble
	X df(n);
	const double fx = function_values[n];
	for (size_type i = 0; i < n; ++i){
		df[i] = (fx- function_values[i]) / h[i];
	}
	return std::make_tuple(fx, df);
}

template <typename F, typename X, typename M>
std::tuple<typename X::value_type, M>
backward_difference::hessian(F f, const X& x)
{
	typedef typename X::value_type value_type;
	typedef typename X::size_type size_type;

	const double eps = std::numeric_limits<value_type>::epsilon();
	auto hmin(std::pow(eps, 1.0/3.0));

	const size_type n = std::distance(x.begin(), x.end());

	X xi(x);
	X xj(x);

	const value_type fx = f(x);
	M H(n, n);

	for (size_type i = 0; i < n; ++i){
		const value_type xii = xi[i];
		const value_type hi = hmin * (1 + fabs(xii));

		// f(x - h_i e_i)
		xi[i] -= hi;
		const value_type fxi = f(xi);
		xi[i] = xii;

		for (size_type j = 0; j <= i; ++j){
			const value_type xjj = xj[j];
			const value_type hj = hmin * (1 + fabs(xjj));

			// f(x - h_j e_j)
			xj[j] -= hj;
			const value_type fxj = f(xj);

			// f(x - h_i e_i - h_j e_j)
			xj[i] -= hi;
			const value_type fxij = f(xj);

			H(i, j) = H(j, i) = (fxij - fxi - fxj + fx)/(hi * hj);

			xj[i] = xii;
			xj[j] = xjj;
		}
	}
	return std::make_tuple(fx, H);
}

} // ns finite_differences
} // ns ook

#endif
