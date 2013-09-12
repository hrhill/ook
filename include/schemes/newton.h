#ifndef OPTIM_SCHEMES_BFGS_H_
#define OPTIM_SCHEMES_BFGS_H_

#include "../more_thuente.h"

namespace ook{

struct newton{

    template <typename F, typename X, typename State>
    static
    State
    initialise(F f, const X& x, const State&)
    {
        return update<F, X, State>(f, x);        
    }

};

} // hs optim

#endif
