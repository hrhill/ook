#define BOOST_TEST_MODULE linear_algebra

#include <iostream>
#include <random>
#include <limits>
#include <ctime>

#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "linear_algebra.hpp"
#include "linear_algebra/special_matrices.hpp"
#include "linear_algebra/norms.hpp"

#include "test_utilities.hpp"

using namespace std;
using namespace linalg;

template <typename Vector, typename Matrix>
int cholesky_tests()
{
    const int n = 5;
    const double threshold = sqrt(std::numeric_limits<double>::epsilon());

    mt19937 rng(std::time(0));
    auto M = generate_spd_matrix<Matrix>(rng, n);

    /// Test solver
    auto x = generate_vector<Vector>(rng, n);
    Vector y(n);
    gemv(1.0, M, x, 0.0, y);

    Vector xsol = cholesky_solve(M, y);

    BOOST_CHECK(norm_infinity(static_cast<const Vector&>(x - xsol)) <= 1e-04);

    /// Test determinant and inversions
    Matrix invM = cholesky_invert(M);
    BOOST_CHECK_CLOSE(cholesky_determinant(M),
                        1.0/cholesky_determinant(invM), threshold);
    BOOST_CHECK_CLOSE(log_cholesky_determinant(M),
                        log(cholesky_determinant(M)), threshold);

    Matrix id1(n, n, 0);
    gemm(1.0, M, invM, 0.0, id1);
    Matrix id2(n, n, 0);
    gemm(1.0, invM, M, 0.0, id2);

    for (int i = 0; i < n; ++i){

        BOOST_CHECK_CLOSE(id1(i, i), 1.0, threshold);
        BOOST_CHECK_CLOSE(id2(i, i), 1.0, threshold);

        for (int j = 0; j < i; ++j){
            BOOST_CHECK(fabs(id1(i, j)) <= threshold);
            BOOST_CHECK(fabs(id1(j, i)) <= threshold);

            BOOST_CHECK(fabs(id2(i, j)) <= threshold);
            BOOST_CHECK(fabs(id2(j, i)) <= threshold);
        }
    }
    // Check the determinants are 1
    BOOST_CHECK_CLOSE(cholesky_determinant(id1), 1.0, threshold);
    BOOST_CHECK_CLOSE(cholesky_determinant(id2), 1.0, threshold);
    return 0;
}

BOOST_AUTO_TEST_CASE(ublas_cholesky_tests)
{
    std::cout << "Testing ublas\n";
    typedef boost::numeric::ublas::vector<double> vector_t;
    typedef boost::numeric::ublas::matrix<double,
            boost::numeric::ublas::column_major> matrix_t;

    int flag = cholesky_tests<vector_t, matrix_t>();
    BOOST_CHECK_EQUAL(flag, 0);
}

#ifdef HAVE_BLAZE
#include <blaze/Math.h>

BOOST_AUTO_TEST_CASE(blaze_cholesky_tests)
{
    std::cout << "Testing blaze\n";
    typedef blaze::DynamicVector<double> vector_t;
    typedef blaze::DynamicMatrix<double, blaze::columnMajor> matrix_t;

    int flag = cholesky_tests<vector_t, matrix_t>();
    BOOST_CHECK_EQUAL(flag, 0);
}

#endif

