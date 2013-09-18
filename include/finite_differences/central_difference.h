#ifndef SMARTODDS_FINITE_DIFFERENCES_CENTRAL_DIFFERENCE_H_
#define SMARTODDS_FINITE_DIFFERENCES_CENTRAL_DIFFERENCE_H_

#include <limits>
#include <cmath>
#include <tuple>

#include <algorithm>
#include "transform.h"

namespace smartodds{
namespace finite_differences{

/// Central difference approximation
struct central_difference{
	/// \brief Calculate the central difference approximation to the gradient of f.
	template <typename F, typename X>
	static 
	std::tuple<typename X::value_type, X>
	gradient(F f, const X& x);

	/// \brief Calculate the central difference approximation to the hessian of f.
	template <typename F, typename X, typename M>
	static 
	std::tuple<typename X::value_type, M>
	hessian(F f, const X& x);
};

template <typename F, typename X>
std::tuple<typename X::value_type, X>
central_difference::gradient(F f, const X& x)
{
	typedef typename X::value_type value_type;
	typedef typename X::size_type size_type;

	// Generate a set of sample points
	// (x + he_1, x + he_2, ..., x+he_n, x-he_1,..., x-he_n, x)
	size_type n = x.size();
	std::vector<X> sample_points(2 * n + 1);
	const value_type hmin(sqrt(std::numeric_limits<value_type>::epsilon()));
	X x_(x);
	X h(n);

	for (size_type i = 0; i < n; ++i)
	{
		const value_type xi = x[i];
		const value_type hx = hmin * (1 + fabs(xi));
		h[i] = hx;

		x_[i] = xi + hx;
		sample_points[i] = x_;
		x_[i] = xi - hx;
		sample_points[i+n] = x_;
		
		x_[i] = xi;
	}
	sample_points[2 * n] = x;
	// evaluate function at each point
	std::vector<value_type> function_values(sample_points.size());
	detail::transform(sample_points.begin(), sample_points.end(), function_values.begin(), f);

	// assemble
	X df(n);
	const value_type fx = function_values[2 * n];	
	for (size_type i = 0; i < n; ++i){
		df[i] = (function_values[i] - function_values[i + n]) / (2.0 * h[i]);
	}
	return std::make_tuple(fx, df);
}

/// Calculate the forward difference approximation to the hessian of f.
template <typename F, typename X, typename M>
std::tuple<typename X::value_type, M>
central_difference::hessian(F f, const X& x)
{
	typedef typename X::value_type value_type;
	typedef typename X::size_type size_type;

	const double eps = std::numeric_limits<value_type>::epsilon();
	auto hmin(std::pow(eps, 1.0/3.0));

	const size_type dim = x.size();
	const value_type fx = f(x);

	X local_x(x);
	M H(dim, dim);
	for (size_type i = 0; i < dim; ++i){
		// Work out diagonal term first

		const value_type hi = hmin * (1.0 + fabs(local_x[i]));
		const value_type xi = local_x[i];

		// f(x + h_i e_i)
		local_x[i] = xi + hi;
		const value_type fp = f(local_x);

		// f(x + 2 h_i e_i)
		local_x[i] = xi + 2.0 * hi;
		const value_type fpp = f(local_x);

		// f(x - h_i e_i)
		local_x[i] = xi - hi;
		const value_type fm = f(local_x);

		// f(x - 2 h_i e_i)
		local_x[i] = xi - 2.0 * hi;
		const value_type fmm = f(local_x);

		H(i, i) = (-fpp + 16.0 * fp - 30.0 * fx + 16.0 * fm - fmm)/(12.0 * hi * hi);

		local_x[i] = xi;

		for (size_type j = 0; j < i; ++j){

			const value_type xj = local_x[j];

			const value_type hj = hmin * (1.0 + fabs(local_x[j]));

			// f(x + h_i e_i + h_j e_j)
			local_x[i] = xi + hi;
			local_x[j] = xj + hj;
			const value_type fpp = f(local_x);

			// f(x - h_i e_i + h_j e_j)
			local_x[i] = xi - hi;
			local_x[j] = xj + hj;
			const value_type fmp = f(local_x);

			// f(x + h_i e_i - h_j e_j)
			local_x[i] = xi + hi;
			local_x[j] = xj - hj;
			const value_type fpm = f(local_x);

			// f(x - h_i e_i - h_j e_j)
			local_x[i] = xi - hi;
			local_x[j] = xj - hj;
			const value_type fmm = f(local_x);

			H(i, j) = H(j, i) = (fpp - fpm - fmp + fmm)/(4.0 * hi * hj);

			local_x[i] = xi;
			local_x[j] = xj;
		}
	}
	return std::make_tuple(fx, H);
}

/// \brief Calculate the central difference approximation to the gradient of f.
/// \details This function computes the central difference difference approximation to the gradient
/// of the scalar function f,
/// \f[ (\nabla f)_i = \frac{f(x + h e_i) - f(x - h e_i)}{2 h} \f].
template <typename F, typename X>
void 
gradient_central_difference(const F& f, const X& x, X& df)
{
	double fx;
	std::tie(fx, df) = central_difference::gradient(f, x);
}

/// \brief User friendly non-reference version.
template <typename F, typename X>
X 
gradient_central_difference(F f, const X& x)
{
	auto dfx = central_difference::gradient(f, x);
	return std::get<1>(dfx);
}

/// \brief Calculate the central difference approximation to the hessian of f.
/// \details This function computes the finite difference approximation to the hessian
/// of the scalar function f,
/// \f{eqnarray*}{ &i \neq j ,(\nabla^2 f)_{ij} & = \frac{f(x + h_x e_i + h_y e_j)
///							      - f(x - h_x e_i + h_y e_j)
///								  - f(x + h_x e_i - h_y e_j)
///								  - f(x - h_x e_i - h_y e_j)}{4 h_i h_j},\\
///  &i = j,  (\nabla^2 f)_{ii} & = \frac{-f(x + 2 h_i e_i)
///												+ 16 f(x + h_i e_i) - 30 f(x)
///												+ 16 f(x - h_i e_i) - f(x - 2 h_i e_i)}{12 h_i^2}
/// \f}
template <typename F, typename X, typename M>
void 
hessian_central_difference(F f, const X& x, M& H)
{
	double fx;
	std::tie(fx, H) = central_difference::hessian<F, X, M>(f, x);
}

} // ns finite_differences
} // ns smartodds

#endif
