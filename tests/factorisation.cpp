/// \file factorisation_test.cpp
#include <iostream>
#include <string>
#include <limits>
#include <random>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/timer.hpp>

#include <boost/mpl/list.hpp>

#define BOOST_TEST_MODULE factorisation_test
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/test_tools.hpp>

#include "factorisations/gmw81.h"

#include "test_utilities.h"

typedef boost::numeric::ublas::matrix<double, boost::numeric::ublas::column_major> matrix_t;

BOOST_AUTO_TEST_CASE(compilation_test){



}
