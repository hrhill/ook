#ifndef SMARTODDS_FINITE_DIFFERENCES_H_
#define SMARTODDS_FINITE_DIFFERENCES_H_

/*! \mainpage finite_differences
 *
 * \section intro_sec Introduction
 *
 * This is a collection of header files containing implementations
 * of basic finite difference calculations. Forward, backward and central difference
 * approximations are available for gradient vector and hessian matrix of scalar
 * valued functions.
 *
 * All approximations can be made to run in parallel by passing the -fopenmp command to the compiler.
 *
 * \section install_sec Installation instructions
 *
 *
 *
 * \section Todo
 * * Intelligent step size calculation.
 * *
 */

#include <smartodds/finite_differences/forward_difference.h>
#include <smartodds/finite_differences/backward_difference.h>
#include <smartodds/finite_differences/central_difference.h>

#endif
