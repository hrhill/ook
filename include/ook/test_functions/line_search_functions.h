#ifndef OOK_TEST_FUNCTIONS_LINE_SEARCH_FUNCTIONS_H_
#define OOK_TEST_FUNCTIONS_LINE_SEARCH_FUNCTIONS_H_

#include <tuple>
#include <cmath>
#include <boost/math/constants/constants.hpp>

namespace ook{

namespace test_functions{

template <typename T>
std::tuple<T, T>
phi51(T a, T b){
    return std::make_tuple(-a/(a * a + b),
                          (a * a - b)/((a * a + b) * (a * a + b)));
}

template <typename T>
std::tuple<T, T>
phi52(T a, T b){
    T t = a + b;
    T t2 = std::pow(t, 2);
    T t3 = std::pow(t, 3);
    T t4 = std::pow(t2, 2);
    T f = t * t4 - 2 * t4;
    T g = 5 * t4 - 8 * t3;
    return std::make_tuple(f, g);
}

template <typename T>
std::tuple<T, T>
phi53(T a, T b, T c){
    const T pi = boost::math::constants::pi<T>();

    T phi0;
    T dphi0;
    if (a <= 1 - b){
        phi0 = 1 - a;
        dphi0 = -1;
    }else if (a >= 1 + b){
        phi0 = a - 1;
        dphi0 = 1;
    }else{
        phi0 = 0.5 * std::pow(a - 1, 2)/b + 0.5 * b;
        dphi0 = (a - 1)/b;
    }
    return std::make_tuple(phi0 + 2 * (1 - b)/(c * pi) * sin(0.5 * c * pi * a),
                           dphi0 + (1 - b) * cos(0.5 * c * pi * a));
}

template <typename T>
std::tuple<T, T>
phi54(T a, T b1, T b2){

    using std::pow;
    T t1 = sqrt(pow(b1, 2) + 1) - b1;
    T t2 = sqrt(pow(b2, 2) + 1) - b2;
    T k1 = sqrt(pow(1 - a, 2) + b2 * b2);
    T k2 = sqrt(pow(a, 2) + pow(b1, 2));
    T f = t1 * k1 + t2 * k2;
    T g = -t1 * (1 - a) / k1 + t2 * a /k2;
    return std::make_tuple(f, g);
}

} // ns test_functions

} // ns ook

#endif
