#define BOOST_TEST_MODULE factorisation_test

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <limits>
#include <random>
#include <string>

#include <boost/mpl/list.hpp>
#include <boost/timer.hpp>

#include "ook.hpp"

using namespace ook;

int
max_magnitude_diagonal_check()
{
    const int n = 5;
    const int max_id = 3;
    const double max_value = 10;

    matrix m = ook::eye(n);
    m(max_id, max_id) = max_value;

    auto x = ook::detail::max_magnitude_diagonal(m);

    BOOST_CHECK_EQUAL(std::get<0>(x), max_id);
    BOOST_CHECK_EQUAL(std::get<1>(x), max_value);
    return 0;
}

int
identity_matrix_check()
{
    const int n = 5;
    matrix A = eye(n);
    ook::detail::gmw81(A);
    for (int i = 0; i < n; ++i)
    {
        BOOST_CHECK_EQUAL(A(i, i), 1.0);
        for (int j = 0; j < i; ++j)
        {
            BOOST_CHECK_EQUAL(A(i, j), 0.0);
        }
    }
    return 0;
}

int
scaled_identity_matrix_check()
{
    const int n = 5;
    std::mt19937 rng(std::time(nullptr));
    auto rnorm = std::bind(std::normal_distribution<>(0, 1), std::ref(rng));
    const double d = exp(rnorm());
    matrix A = d * eye(n);
    ook::detail::gmw81(A);
    for (int i = 0; i < n; ++i)
    {
        BOOST_CHECK_EQUAL(A(i, i), d);
        for (int j = 0; j < i; ++j)
        {
            BOOST_CHECK_EQUAL(A(i, j), 0.0);
        }
    }
    return 0;
}

int
factorisation_equivalence_test()
{
    // Make factorisations agree on positive definite matrices.
    std::mt19937 rng(std::time(nullptr));
    const int n = 5;
    matrix A = sympd(n, rng);
    matrix cholL(A);
    potrf(cholL, 'L');
    matrix gmw81LD(A);
    ook::detail::gmw81(gmw81LD);
    ook::detail::convert_to_cholesky(gmw81LD);

    std::cout << gmw81LD << std::endl;

    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j <= i; ++j)
        {
            BOOST_CHECK_CLOSE(gmw81LD(i, j), cholL(i, j), 1e-01);
        }
    }
    return 0;
}

BOOST_AUTO_TEST_CASE(blaze_factorisation_tests)
{
    max_magnitude_diagonal_check();
    identity_matrix_check();
    scaled_identity_matrix_check();
    factorisation_equivalence_test();
}