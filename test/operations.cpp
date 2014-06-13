#define BOOST_TEST_MODULE operations
#include <iostream>
#include <random>
#include <limits>
#include <ctime>

#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "linear_algebra.hpp"
#include "linear_algebra/operations.hpp"

using namespace std;
using namespace linalg;

template <typename Vector, typename Matrix>
int size_test()
{
    const int n = 5;
    Vector v(n, 0.0);

    BOOST_CHECK_EQUAL(size(v), n);
    return 0;
}

BOOST_AUTO_TEST_CASE(ublas_operations_tests)
{
    std::cout << "Testing ublas\n";
    typedef boost::numeric::ublas::vector<double> vector_t;
    typedef boost::numeric::ublas::matrix<double,
            boost::numeric::ublas::column_major> matrix_t;

    BOOST_CHECK_EQUAL((size_test<vector_t, matrix_t>()), 0);
}


#ifdef HAVE_BLAZE
#include <blaze/Math.h>

BOOST_AUTO_TEST_CASE(blaze_norm_tests)
{
    std::cout << "Testing blaze\n";
    typedef blaze::DynamicVector<double> vector_t;
    typedef blaze::DynamicMatrix<double, blaze::columnMajor> matrix_t;

    BOOST_CHECK_EQUAL((size_test<vector_t, matrix_t>()), 0);
}

#endif
