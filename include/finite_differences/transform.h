#ifndef FINITE_DIFFERENCES_TRANSFORM_H_
#define FINITE_DIFFERENCES_TRANSFORM_H_

#include <algorithm>
#include <parallel/algorithm>

namespace smartodds{
namespace finite_differences{
namespace detail{

template <typename InIter, typename OutIter, typename F>
    void transform(const InIter xbegin, const InIter xend, OutIter ybegin, F f)
{

#ifdef _GLIBCXX_PARALLEL
        __gnu_parallel::_Settings s;
    auto default_strategy = s.algorithm_strategy;
    s.algorithm_strategy = __gnu_parallel::force_parallel;
    __gnu_parallel::_Settings::set(s);
#endif

    std::transform(xbegin, xend, ybegin, f);

#ifdef _GLIBCXX_PARALLEL
    s.algorithm_strategy = default_strategy;
    __gnu_parallel::_Settings::set(s);      
#endif
}

}
} // ns finite differences
} // ns smartodds

#endif