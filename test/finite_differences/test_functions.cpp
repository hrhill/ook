/// \file test_functions.cpp

#include <vector>
#include <unistd.h>
#include <tuple>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

/// Local vector definition.
typedef boost::numeric::ublas::vector<double> vector_t;
typedef boost::numeric::ublas::matrix<double> matrix_t;

/// Rosenbrocks function, minimum is at f(1,1)=0
double 
rosenbrock(const vector_t& x, vector_t& g, matrix_t& H)
{
    sleep(1);
	const double x02 = pow(x(0), 2);
	double f = 100.0 * pow(x(1) - x02, 2) + pow(1.0 - x(0), 2);

	g(0) = -400.0*x(0)*(x(1) - x02) - 2.0*(1.0 - x(0));
	g(1) = 200*(x(1) - x02);

	H(0, 0) = 1200 * x02 -400 * x(1) + 2;
	H(1, 0) = H(0, 1) = -400 * x(0);
	H(1, 1) = 200.0;

    return f;
}

/// Apart from a small region around the optimum this function is flat. Minimum is f(100,100)=-10
double symmetrical_gaussian(const vector_t& x, vector_t& g, matrix_t& H)
{
	const double x0 = x(0) - 100.0;
	const double x1 = x(1) - 100.0;
	const double s = 1500;
	double f = - exp( -pow(x0, 2)/s - pow(x1, 2)/s);

	g(0) = 2 * x0 * f / s;
	g(1) = 2 * x1 * f / s;

	H(0, 0) = - pow(2 * x0 / s, 2) * f + 2/s * f;
	H(1, 0) = H(0, 1) = 4 * x0 * x1 * f / (s * s);
	H(1, 1) = - pow(2 * x1 / s, 2) * f + 2/s * f;

    return f;
}

/// Classical test. Minimum is f(0,0,0,0)=0
double powells_singular(const vector_t& x, vector_t& g, matrix_t& H)
{
    double t1 = ( x(0) + 10.0*x(1) ) * ( x(0) + 10.0 * x(1) );
    double t2 = 5.0*( x(2) - x(3) )*( x(2) - x(3) );
    double t3 = (x(2) - 2.0*x(3))*(x(2) - 2.0*x(3))*(x(2) - 2.0*x(3))*(x(2) - 2.0*x(3));
    double t4 = 10.0*( x(0) - x(3) )*( x(0) - x(3) )*( x(0) - x(3) )*( x(0) - x(3) );
    double f = t1 + t2 + t3 + t4;

    g(0) = 2*(x(0) + 10*x(1)) + 40*(x(0) - x(3))*(x(0) - x(3))*(x(0) - x(3));
    g(1) = 20*(x(0) + 10*x(1));
    g(2) = 10*(x(2) - x(3)) + 4*(x(2) - 2*x(3))*(x(2) - 2*x(3))*(x(2) - 2*x(3));
    g(3) = -10*(x(2) - x(3)) -8*(x(2) - 2*x(3))*(x(2) - 2*x(3))*(x(2) - 2*x(3)) -40*(x(0) - x(3))*(x(0) - x(3))*(x(0) - x(3));

    return f;
}

/// Basic test function.
double paraboloid(const vector_t& x, vector_t& g, matrix_t& H)
{
    g = 2.0 * x;

    for (matrix_t::size_type i = 0; i < H.size1(); ++i){
    	for (matrix_t::size_type j = 0; j < H.size2(); ++j){
    		H(i, j) = 2.0 * (i == j);
    	}
    }
    return pow(boost::numeric::ublas::norm_2(x), 2);
}

double 
rosenbrock_f(const vector_t& x)
{
    sleep(1);
    const double x02 = pow(x(0), 2);
    double f = 100.0 * pow(x(1) - x02, 2) + pow(1.0 - x(0), 2);
    return f;
}

/// Apart from a small region around the optimum this function is flat. Minimum is f(100,100)=-10
double symmetrical_gaussian_f(const vector_t& x)
{
    const double x0 = x(0) - 100.0;
    const double x1 = x(1) - 100.0;
    const double s = 1500;
    double f = - exp( -pow(x0, 2)/s - pow(x1, 2)/s);

    return f;
}

/// Classical test. Minimum is f(0,0,0,0)=0
double powells_singular_f(const vector_t& x)
{
    double t1 = ( x(0) + 10.0*x(1) ) * ( x(0) + 10.0 * x(1) );
    double t2 = 5.0*( x(2) - x(3) )*( x(2) - x(3) );
    double t3 = (x(2) - 2.0*x(3))*(x(2) - 2.0*x(3))*(x(2) - 2.0*x(3))*(x(2) - 2.0*x(3));
    double t4 = 10.0*( x(0) - x(3) )*( x(0) - x(3) )*( x(0) - x(3) )*( x(0) - x(3) );
    double f = t1 + t2 + t3 + t4;

    return f;
}

/// Basic test function.
double paraboloid_f(const vector_t& x)
{
    return pow(boost::numeric::ublas::norm_2(x), 2);
}
