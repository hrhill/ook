#define BOOST_TEST_MODULE finite_differences_test

#include <algorithm>
#include <ctime>
#include <functional>
#include <random>

#include <boost/mpl/list.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/test/test_tools.hpp>
#include <boost/test/unit_test.hpp>

#include "test_functions.hpp"

#include "ook/finite_differences/backward_difference.hpp"
#include "ook/finite_differences/central_difference.hpp"
#include "ook/finite_differences/forward_difference.hpp"

using namespace ook::finite_differences;

template <typename FD, typename F, typename G>
int
checker(F f, G g, int dim)
{
    std::mt19937 rng(std::time(0));
    auto normrnd = bind(std::normal_distribution<>(), std::ref(rng));

    ook::vector x(dim);
    ook::vector df(dim);
    ook::matrix H(dim, dim);
    ook::matrix Hh(dim, dim);

    const int n_tests = 1;

    for (int i = 0; i < n_tests; ++i)
    {
        std::generate(x.begin(), x.end(), normrnd);
        // Calculate finite difference approximation.
        auto dfh = FD::gradient(g, x);
        auto d2fh = FD::template hessian<ook::matrix>(g, x);

        // Evaluate true gradient.
        f(x, df, H);

        // Generate a course upper bound for gradient error
        const double grad_error_bound = 0.1;
        const double hess_error_bound = 1;
        BOOST_REQUIRE_SMALL(ook::norm_inf(static_cast<const ook::matrix&>(
                                H - std::get<1>(d2fh))),
                            hess_error_bound);
        BOOST_REQUIRE_SMALL(ook::norm_inf(static_cast<const ook::vector&>(
                                df - std::get<1>(dfh))),
                            grad_error_bound);
    }
    return 0;
}

typedef boost::mpl::list<forward_difference,
                         backward_difference,
                         central_difference>
    test_types;

BOOST_AUTO_TEST_CASE_TEMPLATE(finite_difference_checker, G, test_types)
{

    BOOST_CHECK_EQUAL(checker<G>(rosenbrock, rosenbrock_f, 2), 0);
    BOOST_CHECK_EQUAL(
        checker<G>(symmetrical_gaussian, symmetrical_gaussian_f, 2), 0);
    BOOST_CHECK_EQUAL(checker<G>(paraboloid, paraboloid_f, 4), 0);
}
