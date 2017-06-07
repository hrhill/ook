#ifndef OOK_FINITE_DIFFERENCE_TEST_FUNCTIONS_HPP_
#define OOK_FINITE_DIFFERENCE_TEST_FUNCTIONS_HPP_

#include "ook/matrix.hpp"
#include "ook/vector.hpp"

double
rosenbrock(const ook::vector& x, ook::vector& g, ook::matrix& H);

double
symmetrical_gaussian(const ook::vector& x, ook::vector& g, ook::matrix& H);

double
powells_singular(const ook::vector& x, ook::vector& g, ook::matrix& H);

double
asymmetrical(const ook::vector& x, ook::vector& g, ook::matrix& H);

double
paraboloid(const ook::vector& x, ook::vector& g, ook::matrix& H);

double
rosenbrock_f(const ook::vector& x);
double
symmetrical_gaussian_f(const ook::vector& x);
double
powells_singular_f(const ook::vector& x);
double
asymmetrical_f(const ook::vector& x);
double
paraboloid_f(const ook::vector& x);

#endif
