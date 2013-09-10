#ifndef OOK_MORE_THUENTE_H_
#define OOK_MORE_THUENTE_H_

#include <tuple>
#include <iomanip>

#include "line_search_conditions.h"
#include "state.h"
#include "dcsrch.h"

namespace ook{

template <typename F>
struct function_evaluator{
  explicit function_evaluator(F function_)
  : function(function_)
  {}

  state
  operator()(state s){
      auto fx = function(s.a);
      ++s.nfev;      
      s.fxap = std::get<0>(fx);
      s.dfxap_dot_p = std::get<1>(fx);
      return s;
  }

  F function;
};

template <typename F>
state
more_thuente_line_search(F phi, state sk, const options& opts){

    dcsrch_struct dcsrch_(sk.fx, sk.dfx_dot_p, sk.a, opts.stpmax - opts.stpmin);

    function_evaluator<F> evaluator(phi);
    sk = evaluator(sk);
    do{
        if (ook::strong_wolfe_conditions(sk, opts.ftol, opts.gtol)){
            sk.value = state_value::convergence;
            break;
        }
        std::tie(sk.value, sk.a) = dcsrch_(sk.a, sk.fxap, sk.dfxap_dot_p, opts);
        sk = evaluator(sk);
    }while (sk.value == state_value::update);

    return sk;
}

} // ns ook

#endif