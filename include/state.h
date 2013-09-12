#ifndef OOK_STATE_H_
#define OOK_STATE_H_

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>

#include "state_value.h"

namespace ook{

struct state{
    typedef double real_type;

    state(const boost::numeric::ublas::vector<double>& x0 = {})
    : 
        x(x0), 
        fx(0.0),
        fxap(0.0),
        dfx_dot_p(0.0),
        dfxap_dot_p(0.0),
        a(0.0),
        nfev(0),        
        dfx(x0.size(), 0.0), 
        dfxap(x0.size(), 0.0), 
        p(x0.size(), 0.0)
    {}

    double fx;
    double fxap;
    double dfx_dot_p;
    double dfxap_dot_p;
    double a;
    state_value value;
    int nfev;
    boost::numeric::ublas::vector<double> x;
    boost::numeric::ublas::vector<double> dfx;
    boost::numeric::ublas::vector<double> dfxap;    
    boost::numeric::ublas::vector<double> p;
};

}  // ns ook

#endif