#ifndef LINALG_TEST_UTILITIES_HPP_
#define LINALG_TEST_UTILITIES_HPP_

#include <random>

#include "linalg/blas/ublas.hpp"
#include "linalg/blas/blaze.hpp"
#include "linalg/operations/ublas.hpp"
#include "linalg/operations/blaze.hpp"

namespace linalg{

template <typename T>
T
generate_matrix(std::mt19937& rng, const int nrow, const int ncol)
{
	T mat(nrow, ncol);
    std::normal_distribution<> n01;
	auto norm_gen = std::bind(n01, std::ref(rng));

	for (int i = 0; i < nrow; ++i){
		for (int j = 0; j < ncol; ++j){
			mat(i, j) = norm_gen();
		}
	}
	return mat;
}

template <typename T>
T
generate_binary_matrix(std::mt19937& rng, const int n, const int m)
{
	T z(n, m, 0.0);
    std::uniform_int_distribution<> uniint(0, m - 1);
	auto idx_gen = std::bind(uniint, std::ref(rng));

	for (int i = 0; i < n; ++i){
		z(i, idx_gen()) = 1.0;
		z(i, idx_gen()) = 1.0;
	}
	return z;
}

template <typename MatrixType, typename VectorType>
MatrixType
diagonal_matrix(const VectorType& d)
{
	auto n = size(d);
	MatrixType m(n, n, 0);
	for (size_t i = 0; i < n; ++i){
		m(i, i) = d(i);
	}
	return m;
}

template <typename MatrixType>
MatrixType
generate_spd_matrix(std::mt19937& rng, const int n)
{
	MatrixType M(n, n);

    std::normal_distribution<> n01;
	auto rnorm = std::bind(n01, std::ref(rng));

	for (int i = 0; i < n; ++i){
		for (int j = 0; j < n; ++j){
			M(i, j) = rnorm();
		}
        M(i, i) = exp(rnorm());
	}
	MatrixType S(n, n);
	gemm(1.0, M, trans(M), 0.0, S);
	return  S;
}

template <typename Vector>
Vector
generate_vector(std::mt19937& rng,
                const int n,
                const double mean = 0,
                const double sd = 1)
{
	Vector v(n);
    std::normal_distribution<> numsig(mean, sd);
	auto rnorm = std::bind(numsig, std::ref(rng));
	std::generate(v.begin(), v.end(), rnorm);
	return v;
}

} // ns linalg

#endif
