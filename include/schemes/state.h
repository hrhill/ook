#ifndef OOK_STATE_H_
#define OOK_STATE_H_

namespace ook{

template <typename RealType, typename VectorType>
struct state{
    enum class status{
        keep_working,
        max_iterations_reached,
        converged
    };

    state()
    {}

    typedef RealType real_type;
    typedef VectorType vector_type;

    real_type fxk;
    real_type alpha;
    real_type fxk_alpha_pk;
    real_type dfxk_dot_pk;    
    real_type dfxk_alpha_pk_dot_pk;    
    vector_type xk;    
    vector_type pk;    
    vector_type dfxk;
    vector_type dfxk_alpha_pk;

    friend
    bool
    operator<(const state& s1, const state& s2)
    {
        return s1.fxk < s2.fxk;
    }

    friend
    bool
    operator>(const state& s1, const state& s2)
    {
        return s2 < s1;
    }

    friend
    bool
    operator<=(const state& s1, const state& s2)
    {
        return !(s1 > s2);
    }

    friend
    bool
    operator>=(const state& s1, const state& s2)
    {
        return !(s1 < s2);
    }
};

}

#endif