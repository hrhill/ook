#ifndef OOK_FINITE_DIFFERENCES_FORWARD_DIFFERENCE_H_
#define OOK_FINITE_DIFFERENCES_FORWARD_DIFFERENCE_H_

#include <limits>
#include <cmath>
#include <tuple>

#include <algorithm>
#include "detail/transform.h"

namespace ook{
namespace finite_differences{

/// Forward difference approximation.
struct forward_difference{
	/// \brief Calculate the forward difference approximation to the gradient of f.
	template <typename F, typename X>
	static
	std::tuple<typename X::value_type, X>
	gradient(F f, const X& x);

	/// \brief Calculate the forward difference approximation to the hessian of f.
	template <typename F, typename X, typename M>
	static
	std::tuple<typename X::value_type, M>
	hessian(F f, const X& x);
};

template <typename F, typename X>
std::tuple<typename X::value_type, X>
forward_difference::gradient(F f, const X& x)
{
	typedef typename X::value_type value_type;
	typedef typename X::size_type size_type;

	// Generate a set of sample points
	// (x + he_1, x + he_2, ..., x)
	size_type n = x.size();
	std::vector<X> sample_points(n + 1);
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
		x_[i] = xi;
	}
	sample_points[n] = x;

	// evaluate function at each point
	std::vector<value_type> function_values(sample_points.size());
	detail::transform(sample_points.begin(), sample_points.end(), function_values.begin(), f);

	// assemble
	X df(n);
	const double fx = function_values[n];
	for (size_type i = 0; i < n; ++i){
	    df[i] = (function_values[i] - fx) / h[i];
	}
	return std::make_tuple(fx, df);
}

template <typename F, typename X, typename M>
std::tuple<typename X::value_type, M>
forward_difference::hessian(F f, const X& x)
{
	typedef typename X::value_type value_type;
	typedef typename X::size_type size_type;

	const value_type eps = std::numeric_limits<value_type>::epsilon();
	value_type hmin(std::pow(eps, 1.0/3.0));

	const size_type dim = x.size();

	X xi(x);
	X xj(x);

	const value_type fx = f(x);
	M H(dim, dim);

#pragma omp parallel for default(none)\
	shared(H, hmin, x)\
	firstprivate(xi, xj, f)

	for (size_type i = 0; i < dim; ++i){
	    const value_type xii = xi[i];
	    const value_type hi = hmin * (1 + fabs(xii));
	    xi[i] += hi;
		// f(x + h_i e_i)
	    const value_type fxi = f(xi);
	    xi[i] = xii;

		for (size_type j = 0; j <= i; ++j){
		    const value_type xjj = xj[j];
		    const value_type hj = hmin * (1 + fabs(xjj));

			// f(x + h_j e_j)
		    xj[j] += hj;
		    const value_type fxj = f(xj);

			// f(x + h_i e_i + h_j e_j)
		    xj[i] += hi;
		    const value_type fxij = f(xj);

#pragma omp critical
		    H(i, j) = H(j, i) = (fxij - fxi - fxj + fx)/(hi * hj);

		    xj[i] = xii;
		    xj[j] = xjj;
		}
	}
	return std::make_tuple(fx, H);
}

/// \brief Calculate the forward difference approximation to the gradient of f.
/// \details This function computes the forward difference difference approximation to the gradient
/// of the scalar function f,
/// \f[ (\nabla f)_i = \frac{f(x + h e_i) - f(x)}{h} \f].
template <typename F, typename X>
void
gradient_forward_difference(F f, const X& x, X& df)
{
	double fx;
	std::tie(fx, df) = forward_difference::gradient(f, x);
}

/// \brief User friendly non-reference version.
template <typename F, typename X>
X
gradient_forward_difference(F f, const X& x)
{
	auto df = forward_difference::gradient(f, x);
	return std::get<1>(df);
}

/// \brief Calculate the forward difference approximation to the hessian of f.
/// \details This function computes the finite difference approximation to the hessian
/// of the scalar function f,
/// \f[ (\nabla^2 f)_{ij} = \frac{f(x + h_i e_i + h_j e_j)
///							      - f(x + h_i e_i)
///								  - f(x + h_j e_j) + f(x)}{h_i h_j} \f].
template <typename F, typename X, typename M>
void 
hessian_forward_difference(F f, const X& x, M& H)
{
	double fx;
	std::tie(fx, H) = forward_difference::hessian<F, X, M>(f, x);
}

/*
 * Check call in client code before un-commenting this, e.g.,
 *
 *  matrix_t H = hessian_forward_difference(f, x);
 *
 *  or
 *  matrix_t H = hessian_forward_difference<matrix_t>(f, x);
 *
template <typename F, typename X, typename M>
auto hessian_forward_difference(F f, const X& x) -> decltype(M()){

	M H(x.size(), x.size(), 0);
	forward_difference fd;
	fd.hessian(f, x, H);
	return H;
}
*/
} // ns finite_differences
} // ns ook

#endif
