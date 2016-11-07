#include <chrono>
#include <iostream>
#include <random>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <boost/numeric/bindings/blas/level1.hpp>
#include <boost/numeric/bindings/blas/level2.hpp>
#include <boost/numeric/bindings/blas/level3.hpp>

#include "ook/newton.h"

namespace ublas = boost::numeric::ublas;

typedef ublas::vector<double> vector_t;
typedef ublas::matrix<double, ublas::column_major> matrix_t;

template <typename F>
void
time_it(F f)
{
    namespace chrono = std::chrono;
    auto t0 = chrono::system_clock::now();
    f();
    auto t1 = chrono::system_clock::now();
    auto diff = chrono::duration_cast<chrono::microseconds>(t1 - t0);
    std::cout << "1 run, taking " << diff.count() / 1e6 << " seconds."
              << std::endl;
}

matrix_t
generate_sdp(std::mt19937& rng, const int n)
{
    auto rnorm = std::bind(std::normal_distribution<>(), rng);
    matrix_t A(n, n);
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            A(i, i) = rnorm();
        }
    }
    matrix_t C(n, n);
    boost::numeric::bindings::blas::gemm(1.0, A, trans(A), 0.0, C);
    return C;
}

vector_t
rnorm(std::mt19937& rng, const int n)
{
    auto rnorm = std::bind(std::normal_distribution<>(), rng);
    vector_t x(n);
    std::generate(x.begin(), x.end(), rnorm);
    return x;
}

template <typename Matrix, typename Vector>
Vector
inline_solve(Matrix A, const Vector& b)
{
    Matrix L = ook::factorisations::gmw81(A);

    // Convert this to a regular cholesky factorised matrix
    // Matrix L = convert_to_cholesky(LD);
    const int n = L.size1();
    for (int j = 0; j < n; ++j)
    {
        const double di = sqrt(L(j, j));
        L(j, j) = di;
        for (int i = j + 1; i < n; ++i)
        {
            L(i, j) *= di;
        }
    }

    boost::numeric::ublas::symmetric_adaptor<Matrix,
                                             boost::numeric::ublas::lower>
        sa(L);
    Matrix b1(b.size(), 1);

    boost::numeric::ublas::column(b1, 0) = b;
    boost::numeric::bindings::lapack::potrs(sa, b1);

    return boost::numeric::ublas::column(b1, 0);
    return b;
}

int
main(int argc, char** argv)
{
    std::mt19937 rng(0);
    int n = atoi(argv[1]);

    auto A = generate_sdp(rng, n);
    auto xsol = rnorm(rng, n);
    vector_t z = ublas::prod(A, xsol);
    /*
        std::cout << "Library\n";
        time_it([A, z, xsol](){
            auto x = ook::detail::solve(A, z);
            std::cout << "error : " << ublas::norm_inf(x - xsol) << std::endl;
        });
    */
    std::cout << std::endl;
    std::cout << "Optimal\n";
    time_it([A, z, xsol]() {
        auto x = inline_solve(A, z);
        std::cout << "error : " << ublas::norm_inf(x - xsol) << std::endl;
    });
}
