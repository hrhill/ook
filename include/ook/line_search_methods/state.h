#ifndef OOK_LINE_SEARCH_METHODS_STATE_H_
#define OOK_LINE_SEARCH_METHODS_STATE_H_

namespace ook{

namespace line_search_method{

namespace detail{

/// \brief State for use with line search method.
template <typename X>
struct 
state{
    typedef X vector_type;
    typedef typename X::value_type value_type;

    state(const int n)
    :
        dfx(n), dfx0(n), p(n), a(1), beta(0)
    {}


    value_type fx;
    vector_type dfx;
    vector_type dfx0;        
    vector_type p;    
    value_type a;    
    value_type beta;      
};

} // ns detail

} // ns line_search_methods

} // ns ook

#endif