#ifndef TRAITS_HPP_
#define TRAITS_HPP_

#include "std_traits.hpp"

namespace linalg{

template <typename T>
using vector_value_type
    = remove_const_reference<decltype(T()[0])>;

template <typename T>
using matrix_value_type
    = remove_const_reference<decltype(T()(0,0))>;

}

#endif
