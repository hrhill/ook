#ifndef OOK_OOK_VECTOR_HPP_
#define OOK_OOK_VECTOR_HPP_

#include <string>
#include <random>

#include <blaze/Math.h>

#include "ook/matrix.hpp"

namespace ook
{
inline namespace v1
{

typedef blaze::DynamicVector<double> vector;

vector n01_vector(int n, std::mt19937& rng);

vector
gaussian_vector(const vector& mu, const matrix& sigma, std::mt19937& rng);

/// \brief \f$ l_1 \f$ norm.
double norm_1(const vector& x);

/// \brief \f$ l_2  \f$ norm.
double norm_2(const vector& x);

/// \brief \f$ l_p  \f$ norm.
double norm_p(const vector& x, int p);

/// \brief \f$ l_{\infty} \f$ norm.
double norm_inf(const vector& x);

/// \brief Read from file.
vector read(int n, const std::string& file);

/// \brief Read from file.
void write(const vector& x, const std::string& file);

} // ns v1
} // ns ook

#endif