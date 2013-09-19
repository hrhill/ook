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

template <typename Matrix>
Matrix
convert_to_L(const Matrix& LD)
{
    int n = LD.size1();
    Matrix L(LD);
    Matrix D(n, n, 0);

    for (int i = 0; i < n; ++i){
        D(i, i) = sqrt(L(i, i));
        L(i, i) = 1.0;
    }
}

BOOST_AUTO_TEST_CASE(compilation_test){

    std::mt19937 rng(std::time(0));
    const int ndim = 3;
    matrix_t A = ook::generate_spd_matrix<matrix_t>(rng, ndim);
    matrix_t cholL = ook::get_lower_cholesky_factor(A);
    matrix_t ldlL = ook::ldlt_factorisation(A);
    matrix_t gmw81L = ook::gmw81(A);    

    std::cout << "chol: " << cholL << std::endl;
    std::cout << "ldl : " << ldlL << std::endl;    
    std::cout << "gmw : " << gmw81L << std::endl;  


}
