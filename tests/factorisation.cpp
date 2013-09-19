/// \file factorisation_test.cpp
#include <iostream>
#include <string>
#include <limits>
#include <random>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>

#include <boost/mpl/list.hpp>

#define BOOST_TEST_MODULE factorisation_test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include "factorisations/cholesky.h"
#include "factorisations/gmw81.h"

#include "test_utilities.h"

typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> matrix_t;

BOOST_AUTO_TEST_CASE(compilation_test){

    std::mt19937 rng(1);
    const int ndim = 3;
    matrix_t A = ook::generate_spd_matrix<matrix_t>(rng, ndim);
    matrix_t cholL = ook::get_lower_cholesky_factor(A);
    matrix_t ldlL = ook::ldlt_factorisation(A);
    matrix_t gmw81L = ook::gmw81(A);    

    std::cout << cholL << std::endl;
    std::cout << ldlL << std::endl;    
    std::cout << gmw81L << std::endl;    
}
