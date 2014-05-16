// Copyright 2013 Harry Hill
//
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

#ifndef OOK_STATE_H_
#define OOK_STATE_H_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>
#include <iomanip>

#include "ook/message.h"
#include "ook/norms.h"
#include "ook/type_traits.h"

namespace ook{
namespace detail{

/// \brief State for use with line search method.
enum class state_tag{init, iterate, final};

template <typename X>
struct
state{
    typedef X vector_type;
    typedef typename X::value_type value_type;
    typedef boost::numeric::ublas::matrix<value_type,
              boost::numeric::ublas::column_major> matrix_type;

    state(const int n = 0, const bool with_matrix = false)
    :
        fx(0),
        dfx_dot_p(0),
        dfx(n),
        dfx0(n),
        p(n),
        dx(n),
        H(n * with_matrix, n * with_matrix, 0.0),
        a(1),
        beta(0),
        iteration(0),
        nfev(0),
        tag(state_tag::init)
    {
        for (int i = 0; i < n * with_matrix; ++i){
            H(i, i) = 1.0;
        }
    }

    friend
    std::ostream&
    operator<<(std::ostream& out, const state& s){
        if (s.tag == state_tag::init){
            out << std::endl
                << std::setw(6) << "n"
                << std::setw(6) << "nfev"
                << std::scientific
                << std::setw(14) << "a"
                << std::setw(14) << "fx"
                << std::setw(14) << "max ||dfx||"
                << std::setw(14) << "max ||dx||" << std::endl;
        }

        if (s.tag == state_tag::iterate){
            out << std::setw(6) << s.iteration
                << std::setw(6) << s.nfev
                << std::scientific
                << std::setw(14) << s.a
                << std::setw(14) << s.fx
                << std::setw(14) << ook::norm_infinity(s.dfx)
                << std::setw(14) << ook::norm_infinity(s.dx);
        }

        if (s.tag == state_tag::final)
        {
            out << "\nstatus : " << s.msg << std::endl;
            out << std::setw(8) << "iter"
                << std::setw(8) << "nfev"
                << std::setw(16) << "fx"
                << std::setw(16) << "max ||dfx||"
                << std::setw(16) << "max ||dx||" << std::endl;

            out << std::setw(8) << s.iteration
                << std::setw(8) << s.nfev
                << std::scientific
                << std::setw(16) << s.fx
                << std::setw(16) << ook::norm_infinity(s.dfx)
                << std::setw(16) << ook::norm_infinity(s.dx) << std::endl;
        }
        return out;
    }

    value_type fx;
    value_type dfx_dot_p;
    vector_type dfx;
    vector_type dfx0;
    vector_type p;
    vector_type dx;
    matrix_type H;
    value_type a;
    value_type beta;
    int iteration;
    int nfev;
    state_tag tag;
    message msg;
};

static_assert(is_regular<state<boost::numeric::ublas::vector<double>>>::value, "state not regular");

} // ns detail

} // ns ook

#endif
