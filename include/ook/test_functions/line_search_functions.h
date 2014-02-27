#ifndef OOK_TEST_FUNCTIONS_LINE_SEARCH_FUNCTIONS_H_
#define OOK_TEST_FUNCTIONS_LINE_SEARCH_FUNCTIONS_H_

#include <tuple>
#include <boost/math/constants/constants.hpp>

namespace ook{

namespace test_functions{

std::tuple<double, double>
phi51(double a, double b){
    return std::make_tuple(-a/(a * a + b),
                          (a * a - b)/((a * a + b) * (a * a + b)));
}

std::tuple<double, double>
phi52(double a, double b){
    double t = a + b;
    double t2 = std::pow(t, 2);
    double t3 = std::pow(t, 3);
    double t4 = std::pow(t2, 2);
    double f = t * t4 - 2 * t4;
    double g = 5 * t4 - 8 * t3;
    return std::make_tuple(f, g);
}

std::tuple<double, double>
phi53(double a, double b, double c){
    const double pi = boost::math::constants::pi<double>();

    double phi0;
    double dphi0;
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

std::tuple<double, double>
phi54(double a, double b1, double b2){

    auto gamma = [](const double b) -> double {
        return sqrt(1.0 + std::pow(b, 2)) - b;
    };

    const double t1 = sqrt(std::pow(1 - a, 2) + pow(b2, 2));
    const double t2 = sqrt(std::pow(a, 2) + std::pow(b1, 2));
    return std::make_tuple(
            gamma(b1) * t1 + gamma(b2) * t2,
          - gamma(b1) * (1 - a)/t1 + gamma(b2) * a/t2);
}

} // ns test_functions

} // ns ook

#endif
