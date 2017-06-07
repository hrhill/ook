// This file is part of ook.
//
// ook is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// ook is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public License
// along with ook.  If not, see <http://www.gnu.org/licenses/>.

#include <iomanip>
#include <iostream>
#include <limits>

#include "ook/message.hpp"
#include "ook/options.hpp"
#include "ook/stream_observer.hpp"

#include "ook/test_functions/more_garbow_hillstrom/rosenbrock.hpp"

#include "ook/bfgs.hpp"
#include "ook/newton.hpp"
#include "ook/nonlinear_cg.hpp"
#include "ook/steepest_descent.hpp"

template <typename F>
struct gradient_only_wrapper
{
    gradient_only_wrapper(F f) : func(f) {}

    std::tuple<double, ook::vector>
    operator()(const ook::vector& x) const
    {
        const int n = x.size();
        double f;
        ook::vector df(n);
        ook::matrix d2f(n, n);

        std::tie(f, df, d2f) = func(x);
        return std::make_tuple(f, df);
    }

    F func;
};

int
main()
{

    ook::options<double> opts;

    typedef ook::test_functions::rosenbrock test_function;
    test_function objective_function;
    gradient_only_wrapper<test_function> wrapper(objective_function);

    ook::vector x(test_function::n, 0.0);
    std::copy(test_function::x0.begin(), test_function::x0.end(), x.begin());
    {
        std::cout << "steepest_descent\n";
        ook::stream_observer<std::ostream> obs(std::cout);
        auto soln = ook::steepest_descent(wrapper, x, opts, obs);
        std::cout << soln.msg << "\n" << soln.x << std::endl;
    }

    {
        std::cout << "fletcher_reeves\n";
        ook::stream_observer<std::ostream> obs(std::cout);
        auto soln = ook::fletcher_reeves(wrapper, x, opts, obs);
        std::cout << soln.msg << "\n" << soln.x << std::endl;
    }

    {
        std::cout << "bfgs\n";
        ook::stream_observer<std::ostream> obs(std::cout);
        auto soln = ook::bfgs(wrapper, x, opts, obs);
        std::cout << soln.msg << "\n" << soln.x << std::endl;
    }

    {
        std::cout << "newton\n";
        ook::stream_observer<std::ostream> obs(std::cout);
        auto soln = ook::newton(objective_function, x, opts, obs);
        std::cout << soln.msg << "\n" << soln.x << std::endl;
    }
}
