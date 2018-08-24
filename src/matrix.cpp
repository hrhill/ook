#include "ook/matrix.hpp"
#include <functional>

namespace ook
{
inline namespace v1
{
matrix
eye(const size_t n)
{
    matrix id(n, n, 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        id(i, i) = 1.0;
    }
    return id;
}

matrix
n01_matrix(size_t m, size_t n, std::mt19937& rng)
{
    auto rnorm = std::bind(std::normal_distribution<>(), std::ref(rng));
    matrix a(m, n, 0.0);
    for (size_t i = 0; i < m; ++i)
    {
        for (size_t j = 0; j < n; ++j)
        {
            a(i, j) = rnorm();
        }
    }
    return a;
}

matrix
sympd(size_t n, std::mt19937& rng)
{
    auto rnorm = std::bind(std::normal_distribution<>(), std::ref(rng));
    matrix a(n, n, 0.0);
    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < i; ++j)
        {
            a(i, j) = rnorm();
        }
        a(i, i) = std::pow(rnorm(), 2);
    }
    return a * trans(a);
}

double
norm_inf(const matrix& a)
{
    double nx = 0;
    for (size_t i = 0; i < a.rows(); ++i)
    {
        double rs = 0.0;
        for (size_t j = 0; j < a.columns(); ++j)
        {
            rs += a(i, j);
        }
        if (rs > nx)
            nx = rs;
    }
    return nx;
}

double
norm_inf(const symmetric_matrix& a)
{
    double nx = 0;
    for (size_t i = 0; i < a.rows(); ++i)
    {
        double rs = 0.0;
        for (size_t j = 0; j < a.columns(); ++j)
        {
            rs += a(i, j);
        }
        if (rs > nx)
            nx = rs;
    }
    return nx;
}
}
}
