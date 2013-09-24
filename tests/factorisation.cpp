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

#include "ook/factorisations/cholesky.h"
#include "ook/factorisations/ldlt.h"
#include "ook/factorisations/gmw81.h"

#include "test_utilities.h"

typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> matrix_t;

matrix_t
convert_to_cholesky(const matrix_t& LD)
{
    int n = LD.size1();
    matrix_t L(LD);
    matrix_t D(n, n, 0);

    for (int i = 0; i < n; ++i){
        D(i, i) = sqrt(L(i, i));
        L(i, i) = 1.0;
    }
    return boost::numeric::ublas::prod(L, D);
}

BOOST_AUTO_TEST_CASE(identity_matrix_check){

    const int ndim = 5;

    matrix_t A = boost::numeric::ublas::identity_matrix<double>(ndim);
    matrix_t L = ook::factorisations::gmw81(A);
    const int n = L.size1();
    for (int i = 0; i < n; ++i){
        BOOST_CHECK_EQUAL(L(i, i), 1.0);
        for (int j = 0; j < i; ++j){
            BOOST_CHECK_EQUAL(L(i, j), 0.0);
        }
    }
}

BOOST_AUTO_TEST_CASE(scaled_identity_matrix_check){

    const int ndim = 5;
    std::mt19937 rng(std::time(0));
    auto rnorm = std::bind(std::normal_distribution<>(0, 1), std::ref(rng));
    const double d = exp(rnorm());

    matrix_t A = d * boost::numeric::ublas::identity_matrix<double>(ndim);
    matrix_t L = ook::factorisations::gmw81(A);
    const int n = L.size1();
    for (int i = 0; i < n; ++i){
        BOOST_CHECK_EQUAL(L(i, i), d);
        for (int j = 0; j < i; ++j){
            BOOST_CHECK_EQUAL(L(i, j), 0.0);
        }
    }
}

BOOST_AUTO_TEST_CASE(factorisation_equivalence_test){

    // Make all three factorisations agree on positive definite matrices.
    std::mt19937 rng(std::time(0));
    const int ndim = 5;
    matrix_t A = ook::generate_spd_matrix<matrix_t>(rng, ndim);

    matrix_t cholL = ook::factorisations::cholesky(A);
    matrix_t ldlLD = ook::factorisations::ldlt(A);
    matrix_t gmw81LD = ook::factorisations::gmw81(A);    

    matrix_t ldlL = convert_to_cholesky(ldlLD);
    matrix_t gmw81L = convert_to_cholesky(gmw81LD);    

    for (int i = 0; i < ndim; ++i){
        for (int j = 0; j <= i; ++j){
            BOOST_CHECK_CLOSE(ldlLD(i, j), gmw81LD(i, j), 1e-01);            
            BOOST_CHECK_CLOSE(ldlL(i, j), gmw81L(i, j), 1e-01);
            BOOST_CHECK_CLOSE(gmw81L(i, j),cholL(i, j), 1e-01);            
        }
    }
}
