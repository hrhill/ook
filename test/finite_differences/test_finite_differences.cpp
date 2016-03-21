#define BOOST_TEST_MODULE finite_differences_test

#include <algorithm>
#include <ctime>
#include <random>
#include <functional>

#include <boost/test/unit_test.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/bind.hpp>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "test_functions.hpp"

#include "ook/finite_differences/forward_difference.hpp"
#include "ook/finite_differences/backward_difference.hpp"
#include "ook/finite_differences/central_difference.hpp"

using namespace ook::finite_differences;

typedef boost::numeric::ublas::vector<double> vector_t;
typedef boost::numeric::ublas::matrix<double> matrix_t;


template <typename FD, typename F, typename G>
int checker(F f, G g, int dim){

    namespace ublas = boost::numeric::ublas;

    std::mt19937 rng(std::time(0));
    auto normrnd = bind(std::normal_distribution<>(), std::ref(rng));

    vector_t x(dim);
    vector_t df(dim);
    matrix_t H(dim, dim);
    matrix_t Hh(dim, dim);

    const int n_tests = 1;

    for (int i = 0; i < n_tests; ++i){
        std::generate(x.begin(), x.end(), normrnd);
        // Calculate finite difference approximation.
        auto dfh = FD::gradient(g, x);
        auto d2fh = FD::template hessian<matrix_t>(g, x);

        // Evaluate true gradient.
        f(x, df, H);

        // Generate a course upper bound for gradient error
        const double grad_error_bound = 0.1;
        const double hess_error_bound = 1;
        BOOST_REQUIRE_SMALL(static_cast<double>(ublas::norm_inf(H - std::get<1>(d2fh))), hess_error_bound);
        BOOST_REQUIRE_SMALL(ublas::norm_inf(df - std::get<1>(dfh)), grad_error_bound);
    }
    return 0;
}

typedef boost::mpl::list<forward_difference,
                          backward_difference,
                          central_difference
                          > test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(finite_difference_checker, G, test_types){

    BOOST_CHECK_EQUAL(checker<G>(rosenbrock, rosenbrock_f, 2), 0);
    BOOST_CHECK_EQUAL(checker<G>(symmetrical_gaussian, symmetrical_gaussian_f, 2), 0);
    BOOST_CHECK_EQUAL(checker<G>(paraboloid, paraboloid_f, 4), 0);
}

