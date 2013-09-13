#ifndef OOK_STATE_H_
#define OOK_STATE_H_

#include <iostream>
#include <boost/numeric/ublas/vector.hpp>

#include "state_value.h"

namespace ook{

template <typename T>
struct state_traits{
    typedef typename T::value_type real_type;
    
    inline
    static 
    typename T::size_type
    size(const T& t){
        return t.size();
    }
};

template <>
struct state_traits<double>{
    typedef double real_type;
    
    // zero initialize a double
    static
    size_t
    size(const double& d){
        return size_t(0);
    }
};

template <typename T>
struct state{
    typedef typename state_traits<T>::real_type real_type;

    explicit state(const T& x0 = T())
    : 
        fx(0.0),
        dfx_dot_p(0.0),
        a(0.0),
        value(state_value::start),        
        x(x0), 
        dfx(state_traits<T>::size(x0))
    {}

    real_type fx;
    real_type dfx_dot_p;
    real_type a;
    state_value value;
    T x;
    T dfx;

    friend
    std::ostream&
    operator<<(std::ostream& out, const state& s){
        return out << "{ fx : " << s.fx << ","
                   << " a : " << s.a << ", "
                   << " dfx_dot_p : " << s.dfx_dot_p << ", "
                   << " value : " << s.value << "} ";
    }    
};

}  // ns ook

#endif