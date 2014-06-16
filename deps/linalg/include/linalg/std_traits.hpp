#ifndef REMOVE_CONSTANT_REFERENCE_HPP_
#define REMOVE_CONSTANT_REFERENCE_HPP_

#include <type_traits>

template <typename T>
using remove_const_reference
    = typename std::remove_const<
        typename std::remove_reference<T>::type>::type;

# endif
