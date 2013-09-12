#ifndef OOK_STATE_H_
#define OOK_STATE_H_

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>

#include "state_value.h"

namespace ook{

struct state{
    typedef double real_type;

    state()
    : 
        nfev(0)
    {}

    double fx;
    double fxap;
    double dfx_dot_p;
    double dfxap_dot_p;
    double a;
    state_value value;
    int nfev;
    boost::numeric::ublas::vector<double> x;
    boost::numeric::ublas::vector<double> p;
};

}  // ns ook

#endif