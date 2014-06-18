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

#ifndef OOK_TYPE_TRAITS_H_
#define OOK_TYPE_TRAITS_H_

#include <type_traits>

namespace ook{

template<typename T>
struct is_regular
: std::integral_constant<bool,
    std::is_default_constructible<T>::value &&
    std::is_copy_constructible<T>::value &&
    std::is_move_constructible<T>::value &&
    std::is_copy_assignable<T>::value &&
    std::is_move_assignable<T>::value >
{};

}

#endif
