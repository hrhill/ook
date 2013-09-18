/// \file utilities.h
#ifndef SMARTODDS_FILTERING_TEMPLATE_LIBRARY_TEST_UTILITIES_H_
#define SMARTODDS_FILTERING_TEMPLATE_LIBRARY_TEST_UTILITIES_H_

#include <random>

#include <smartodds/finite_differences.h>
#include "ftl/linear_algebra.h"
#include "ftl/types.h"

namespace ftl{

template <typename T>
T
generate_matrix(std::mt19937& rng, const int nrow, const int ncol){

	T mat(nrow, ncol);
	auto norm_gen = std::bind(std::normal_distribution<>(), std::ref(rng));

	for (int i = 0; i < nrow; ++i){
		for (int j = 0; j < ncol; ++j){
			mat(i, j) = norm_gen();
		}
	}
	return mat;
}

template <typename T>
T
generate_binary_matrix(std::mt19937& rng, const int n, const int m){

	T z(n, m, 0.0);

	auto idx_gen = std::bind(std::uniform_int_distribution<>(0.0, m - 1), std::ref(rng));

	for (int i = 0; i < n; ++i){
		z(i, idx_gen()) = 1.0;
		z(i, idx_gen()) = 1.0;
	}
	return z;
}

template <typename MatrixType>
MatrixType
constant_diagonal_matrix(const int n, const double d){
	MatrixType m(n, n, 0.0);
	for (int i = 0; i < n; ++i){
		m(i, i) = d;
	}
	return m;
}

template <typename MatrixType>
MatrixType
identity_matrix(const int n){
	return constant_diagonal_matrix<MatrixType>(n, 1.0);
}

template <typename MatrixType, typename VectorType>
MatrixType
diagonal_matrix(const VectorType& d){

	auto n = d.size();
	MatrixType m(n, n, 0);
	for (typename MatrixType::size_type i = 0; i < n; ++i){
		m(i, i) = d(i);
	}
	return m;
}

template <typename MatrixType>
MatrixType
generate_spd_matrix(std::mt19937& rng, const int n){

	MatrixType M(n, n, 0.0);

	auto norm_gen = std::bind(std::normal_distribution<>(), std::ref(rng));

	for (int i = 0; i < n; ++i){
		for (int j = 0; j <= i; ++j){
			M(i, j) = norm_gen();
		}
	}
	MatrixType S(n, n);
	boost::numeric::bindings::blas::gemm(1.0, M, trans(M), 0.0, S);
	return  S;
}

template <typename Vector>
Vector
generate_vector(std::mt19937& rng, const int n, const double mean = 0, const double sd = 1)
{
	Vector v(n);
	auto norm_gen = std::bind(std::normal_distribution<>(mean, sd), std::ref(rng));
	std::generate(v.begin(), v.end(), norm_gen);
	return v;
}

} // ns ftl

#endif
