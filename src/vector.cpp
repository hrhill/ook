#include "ook/vector.hpp"
#include "ook/matrix.hpp"

#include <fstream>

namespace ook
{
inline namespace v1
{
vector
n01_vector(int n, std::mt19937& rng)
{
    auto rnorm = std::bind(std::normal_distribution<>(), std::ref(rng));
    vector x(n, 0);
    std::generate(x.begin(), x.end(), rnorm);
    return x;
}

vector
gaussian_vector(const vector& mu, const matrix& sigma, std::mt19937& rng)
{
    auto x = n01_vector(mu.size(), rng);
    matrix L;
    llh(sigma, L); // LLH decomposition
    return mu + L * x;
}

/// \brief \f$ l_1 \f$ norm.
double
norm_1(const vector& x)
{
    double nx = 0.0;
    for (const auto& xi : x) {
        nx += fabs(xi);
}
    return nx;
}

/// \brief \f$ l_2  \f$ norm.
double
norm_2(const vector& x)
{
    double nx = 0.0;
    for (const auto& xi : x) {
        nx += xi * xi;
}
    return sqrt(nx);
}

/// \brief \f$ l_p  \f$ norm.
double
norm_p(const vector& x, int p)
{
    assert(p > 0);
    double nx = 0.0;
    for (const auto& xi : x) {
        nx += std::pow(xi, p);
}
    return exp(log(nx) / p);
}

/// \brief \f$ l_{\infty} \f$ norm.
double
norm_inf(const vector& x)
{
    double nx(0.0);
    for (const auto& xi : x)
    {
        const double fxi = fabs(xi);
        if (fxi > nx) {
            nx = fxi;
}
    }
    return nx;
}

/// \brief Read from file.
vector
read(int n, const std::string& file)
{
    std::ifstream in(file);
    if (!in.good())
    {
        throw std::runtime_error("Cannot open " + file);
    }
    vector x(n, 0.0);
    for (int i = 0; i < n; ++i)
    {
        in >> x[i];
    }
    return x;
}

/// \brief Write to file.
void
write(const vector& x, const std::string& file)
{
    std::ofstream out(file);
    for (const auto& xi : x)
    {
        out << xi << " ";
    }
}
}  // namespace v1
}  // namespace ook
