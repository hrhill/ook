#include <iostream>
#include <iomanip>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/io.hpp>

#include "ook/factorisations/gmw81.h"

typedef boost::numeric::ublas::matrix<double> matrix_t;
typedef boost::numeric::ublas::vector<double> vector_t;

void pretty_print(matrix_t m){

    for (matrix_t::size_type i = 0; i < m.size1(); ++i){
        std::cout << "|";
        for (matrix_t::size_type j = 0; j < m.size2(); ++j){
            std::cout << std::setprecision(4) << std::setw(8) << m(i, j);
        }
        std::cout << "|" << std::endl;
    }
    std::cout << std::endl;
}

int main(){

    std::cout << 
        "\nThis is code to check that the implementation of gmw81\n"
        "gives the same numbers as the example on page 111 of\n"
        "\"Practical Optimization\" by Gill, Murray and Wright, which takes\n"
        "the matrix from example 4.7 (p109) and applies gmw81 to it. The\n"
        "input matrix G is given by,\n\n"
        "|       1       1       2|\n"
        "|       1  1 + 10e-20   3|\n"
        "|       2       3       1|\n";

    const int n = 3;
    matrix_t G(n, n);
    matrix_t L(n, n);    
    matrix_t D(n, n);

    G <<= 1.0, 1.0,         2.0, 
          1.0, nextafter(1.0, 2), 3.0,
          2.0, 3.0,         1.0;

    std::cout << "\nThe result of the algorithm is a lower triangular matrix L and a\n"
                 "diagonal matrix D, given by,\n\n";
    L <<= 1.0,    0.0,    0.0,
          0.2652, 1.0,    0.0,
          0.5304, 0.4295, 1.0;
    D <<= 3.771,   0.0, 0.0,
          0.0,   5.750, 0.0,
          0.0,     0.0, 1.121;

    std::cout << "L = \n";
    pretty_print(L);
    std::cout << "D = \n";    
    pretty_print(D);    

    std::cout << "Output from gmw81\n";
    matrix_t LD = ook::factorisations::gmw81(G);
    matrix_t d(n, n, 0.0);
    for (int i = 0; i < n; ++i){
        d(i, i) = LD(i, i);
        LD(i, i) = 1.0;
    }
    std::cout << "L = \n";    
    pretty_print(LD);
    std::cout << "D = \n";    
    pretty_print(d);
}

