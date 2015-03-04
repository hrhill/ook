/// \file test_func.h

#ifndef OOK_FINITE_DIFFERENCE_TEST_FUNCTIONS_H
#define OOK_FINITE_DIFFERENCE_TEST_FUNCTIONS_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>


double rosenbrock(const boost::numeric::ublas::vector<double>& x,
						boost::numeric::ublas::vector<double>& g,
						boost::numeric::ublas::matrix<double>& H);

double symmetrical_gaussian(const boost::numeric::ublas::vector<double>& x,
								  boost::numeric::ublas::vector<double>& g,
								  boost::numeric::ublas::matrix<double>& H);

double powells_singular(const boost::numeric::ublas::vector<double>& x,
							  boost::numeric::ublas::vector<double>& g,
							  boost::numeric::ublas::matrix<double>& H);

double asymmetrical(const boost::numeric::ublas::vector<double>& x,
						  boost::numeric::ublas::vector<double>& g,
						  boost::numeric::ublas::matrix<double>& H);

double paraboloid(const boost::numeric::ublas::vector<double>& x,
						boost::numeric::ublas::vector<double>& g,
						boost::numeric::ublas::matrix<double>& H);

double rosenbrock_f(const boost::numeric::ublas::vector<double>& x);
double symmetrical_gaussian_f(const boost::numeric::ublas::vector<double>& x);
double powells_singular_f(const boost::numeric::ublas::vector<double>& x);
double asymmetrical_f(const boost::numeric::ublas::vector<double>& x);
double paraboloid_f(const boost::numeric::ublas::vector<double>& x);

#endif
