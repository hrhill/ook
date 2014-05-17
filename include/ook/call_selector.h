// Copyright 2013 Harry Hill
//
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

#ifndef OOK_CALL_SELECTOR_H_
#define OOK_CALL_SELECTOR_H_

#include <tuple>

namespace ook{
namespace detail{

/// \brief Meta function to select the right function call
/// based on the properties of the return type.
template <int dim>
struct call_selector{};

template <>
struct call_selector<1>
{
    template <typename F, typename X, typename State>
    static
    void
    call(F f, const X& x, State& s){
        ++s.nfev;
        std::tie(s.fx) = f(x);
    }
};

template <>
struct call_selector<2>
{
    template <typename F, typename X, typename State>
    static
    void
    call(F f, const X& x, State& s){
        ++s.nfev;
        std::tie(s.fx, s.dfx) = f(x);
    }
};

template <>
struct call_selector<3>{
    template <typename F, typename X, typename State>
    static
    void
    call(F f, const X& x, State& s){
        ++s.nfev;
        std::tie(s.fx, s.dfx, s.H) =  f(x);
    }
};

} // ns detail
} // ns ook

#endif
