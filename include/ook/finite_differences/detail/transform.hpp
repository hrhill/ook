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

#ifndef OOK_FINITE_DIFFERENCES_DETAIL_TRANSFORM_HPP_
#define OOK_FINITE_DIFFERENCES_DETAIL_TRANSFORM_HPP_

#include <algorithm>
#include <boost/config.hpp>

#if BOOST_GCC
#include <parallel/algorithm>
#endif

namespace ook{
namespace finite_differences{
namespace detail{

template <typename InIter, typename OutIter, typename F>
void transform(const InIter xbegin, const InIter xend, OutIter ybegin, F f)
{
/*
#ifdef _GLIBCXX_PARALLEL
        __gnu_parallel::_Settings s;
    auto default_strategy = s.algorithm_strategy;
    s.algorithm_strategy = __gnu_parallel::force_parallel;
    __gnu_parallel::_Settings::set(s);
#endif
*/
    std::transform(xbegin, xend, ybegin, f);
/*
#ifdef _GLIBCXX_PARALLEL
    s.algorithm_strategy = default_strategy;
    __gnu_parallel::_Settings::set(s);
#endif
*/
}

}
} // ns finite differences
} // ns ook

#endif
