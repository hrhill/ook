#ifndef OOK_LINE_SEARCH_METHODS_STATE_H_
#define OOK_LINE_SEARCH_METHODS_STATE_H_

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <iostream>
#include <iomanip>

#include "ook/message.h"

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

    state(const int n, const bool with_matrix = false)
    :
         dfx(n),
         dfx0(n),
         p(n),
         dx(n),
         H(n * with_matrix, n * with_matrix, 0.0),
         a(1),
         beta(0),
         iteration(0),
         tag(state_tag::init)
    {
        if (with_matrix){
            for (int i = 0; i < n; ++i){
                H(i, i) = 1.0;
            }
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
                << std::setw(6) << 0//s.nfev_total
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
                << std::setw(8) << 0//s.nfev_total
                << std::scientific
                << std::setw(16) << s.fx
                << std::setw(16) << ook::norm_infinity(s.dfx)
                << std::setw(16) << ook::norm_infinity(s.dx) << std::endl;
        }
        return out;
    }

    value_type fx;
    vector_type dfx;
    vector_type dfx0;
    vector_type p;
    vector_type dx;
    matrix_type H;
    value_type a;
    value_type beta;
    int iteration;
    state_tag tag;
    message msg;
};

} // ns detail

} // ns ook

#endif
