#include <gtest/gtest.h>

#include <gmpxx.h>

#include "quadratic.hpp"

TEST(QuadraticSieve, FactorsSmallSemiprime)
{
    const mpz_class p = 10007;
    const mpz_class q = 10009;
    const mpz_class N = p * q;

    const auto factors = primetools::QuadraticSieveFactor(N);
    ASSERT_TRUE(factors.has_value());
    EXPECT_NE(factors->first, 1);
    EXPECT_NE(factors->second, 1);
    EXPECT_EQ(factors->first * factors->second, N);
    EXPECT_TRUE(factors->first == p || factors->first == q);
    EXPECT_TRUE(factors->second == p || factors->second == q);
}

TEST(QuadraticSieve, ReturnsNulloptOnPrime)
{
    const mpz_class N = 10007;
    const auto factors = primetools::QuadraticSieveFactor(N);
    EXPECT_FALSE(factors.has_value());
}

TEST(QuadraticSieve, FactorsLargerSemiprime)
{
    // Two well-known 32-bit primes.
    const mpz_class p = 1000000007;
    const mpz_class q = 1000000009;
    const mpz_class N = p * q;

    const auto factors = primetools::QuadraticSieveFactor(N);
    ASSERT_TRUE(factors.has_value());
    EXPECT_EQ(factors->first * factors->second, N);
    EXPECT_TRUE(factors->first == p || factors->first == q);
    EXPECT_TRUE(factors->second == p || factors->second == q);
}

TEST(QuadraticSieve, HandlesEvenInput)
{
    const mpz_class N = mpz_class(2) * mpz_class(263);
    const auto factors = primetools::QuadraticSieveFactor(N);
    ASSERT_TRUE(factors.has_value());
    EXPECT_EQ(factors->first * factors->second, N);
    EXPECT_TRUE(factors->first == 2 || factors->second == 2);
}

TEST(QuadraticSieve, FactorsPerfectSquare)
{
    const mpz_class a = 12345;
    const mpz_class N = a * a;
    const auto factors = primetools::QuadraticSieveFactor(N);
    ASSERT_TRUE(factors.has_value());
    EXPECT_EQ(factors->first, a);
    EXPECT_EQ(factors->second, a);
}

TEST(QuadraticSieve, RecommendedBIsSane)
{
    const mpz_class smallN = mpz_class(10007) * mpz_class(10009);
    const mpz_class largerN = (mpz_class(1) << 140) + 12345;

    const uint32_t b1 = primetools::QuadraticSieveRecommendedFactorBaseBound(smallN);
    const uint32_t b2 = primetools::QuadraticSieveRecommendedFactorBaseBound(largerN);

    EXPECT_GT(b1, 0u);
    EXPECT_GT(b2, 0u);
    EXPECT_GE(b2, b1);
}
