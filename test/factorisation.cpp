/// \file factorisation_test.cpp
#include <iostream>
#include <string>
#include <limits>
#include <random>

#include <boost/timer.hpp>

#include <boost/mpl/list.hpp>

#define BOOST_TEST_MODULE factorisation_test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include "linear_algebra/blas.hpp"
#include "linear_algebra/operations.hpp"
#include "linear_algebra/factorisations/cholesky.hpp"
#include "linear_algebra/factorisations/ldlt.hpp"
#include "linear_algebra/factorisations/gmw81.hpp"
#include "linear_algebra/special_matrices.hpp"

#include "test_utilities.hpp"

template <typename Matrix>
Matrix
convert_to_cholesky(Matrix LD)
{
    int n = LD.size1();
    Matrix L(LD);
    Matrix D(n, n, 0);

    for (int i = 0; i < n; ++i){
        D(i, i) = sqrt(L(i, i));
        L(i, i) = 1.0;
    }
    return boost::numeric::ublas::prod(L, D);
}

template <typename Matrix>
int identity_matrix_check()
{
    const int ndim = 5;

    Matrix A = linalg::identity_matrix<Matrix>(ndim);
    Matrix L = linalg::factorisations::gmw81(A);
    const int n = linalg::num_rows(L);
    for (int i = 0; i < n; ++i){
        BOOST_CHECK_EQUAL(L(i, i), 1.0);
        for (int j = 0; j < i; ++j){
            BOOST_CHECK_EQUAL(L(i, j), 0.0);
        }
    }
    return 0;
}

template <typename Matrix>
int scaled_identity_matrix_check()
{
    const int ndim = 5;
    std::mt19937 rng(std::time(0));
    auto rnorm = std::bind(std::normal_distribution<>(0, 1), std::ref(rng));
    const double d = exp(rnorm());

    Matrix A = d * linalg::identity_matrix<Matrix>(ndim);
    Matrix L = linalg::factorisations::gmw81(A);
    const int n = linalg::num_rows(L);
    for (int i = 0; i < n; ++i){
        BOOST_CHECK_EQUAL(L(i, i), d);
        for (int j = 0; j < i; ++j){
            BOOST_CHECK_EQUAL(L(i, j), 0.0);
        }
    }
    return 0;
}

template <typename Matrix>
int factorisation_equivalence_test()
{
    // Make all three factorisations agree on positive definite matrices.
    std::mt19937 rng(std::time(0));
    const int ndim = 5;
    Matrix A = linalg::generate_spd_matrix<Matrix>(rng, ndim);

    Matrix cholL = linalg::factorisations::cholesky(A);
    Matrix ldlLD = linalg::factorisations::ldlt(A);
    Matrix gmw81LD = linalg::factorisations::gmw81(A);

    Matrix ldlL = convert_to_cholesky(ldlLD);
    Matrix gmw81L = convert_to_cholesky(gmw81LD);

    for (int i = 0; i < ndim; ++i){
        for (int j = 0; j <= i; ++j){
            BOOST_CHECK_CLOSE(ldlLD(i, j), gmw81LD(i, j), 1e-01);
            BOOST_CHECK_CLOSE(ldlL(i, j), gmw81L(i, j), 1e-01);
            BOOST_CHECK_CLOSE(gmw81L(i, j),cholL(i, j), 1e-01);
        }
    }
    return 0;
}


BOOST_AUTO_TEST_CASE(ublas_factorisation_tests)
{
    std::cout << "Testing ublas\n";
    typedef boost::numeric::ublas::matrix<double,
            boost::numeric::ublas::column_major> matrix_t;

    int flag = factorisation_equivalence_test<matrix_t>();
    BOOST_CHECK_EQUAL(flag, 0);
}
/*
#ifdef HAVE_BLAZE
#include <blaze/Math.h>

BOOST_AUTO_TEST_CASE(blaze_factorisation_tests)
{
    std::cout << "Testing blaze\n";
    typedef blaze::DynamicMatrix<double, blaze::columnMajor> matrix_t;

    int flag = factorisation_equivalence_test<matrix_t>();
    BOOST_CHECK_EQUAL(flag, 0);
}

#endif
*/
