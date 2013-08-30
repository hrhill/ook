/// \file optimise_test.cpp
#include <iostream>
#include <string>
#include <limits>
#include <random>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>

#include <boost/mpl/list.hpp>

#define BOOST_TEST_MODULE optimise
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include "optimise.h"
#include "schemes/newton.h"

typedef boost::mpl::list<
                ook::newton//, optim::bfgs, optim::lbfgs
                > optimiser_types;

std::tuple<double, boost::numeric::ublas::vector<double>>
objective_function(const boost::numeric::ublas::vector<double>& x)
{
    return std::make_tuple(0.5 * std::pow(boost::numeric::ublas::norm_2(x), 2), x);
}

BOOST_AUTO_TEST_CASE_TEMPLATE(optimiser_compilation_check, T, optimiser_types){

    boost::numeric::ublas::vector<double> x0(5, 1.1);
    typename T::options_type opts;
    typedef ook::state<double, boost::numeric::ublas::vector<double>> state_type;
    auto soln = ook::optimise<T, state_type>(objective_function, x0, opts);
}

BOOST_AUTO_TEST_CASE(compilation_test){
    BOOST_CHECK(false);
}
