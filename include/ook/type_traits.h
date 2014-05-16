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
