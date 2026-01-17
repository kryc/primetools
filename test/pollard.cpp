#include <array>
#include <limits>
#include <vector>

#include <gtest/gtest.h>

#include "pollard.hpp"

// 32-bit primes that are close together
std::vector<std::array<uint32_t, 2>> PollardTestCases = {
    { { 31699, 65521 } },
    { { 919, 60029 } },
    { { 2473, 5581 } },
    { { 248789, 208319 } }
};

TEST(Pollard, BrentPollardsRho)
{
    for (const auto& testCase : PollardTestCases) {
        const mpz_class p = testCase[0];
        const mpz_class q = testCase[1];
        const mpz_class N = p * q;

        const auto factors = primetools::BrentPollardsRho<mpz_class>(N);
        ASSERT_TRUE(factors.has_value());
        EXPECT_TRUE(factors->first == p || factors->first == q);
        EXPECT_TRUE(factors->second == p || factors->second == q);
    }
}

TEST(Pollard, PollardsPMinus1)
{
    // Multiple deterministic semiprimes where p-1 is B-smooth.
    struct Case {
        uint64_t p;
        uint64_t q;
        size_t B;
        size_t bases;
    };

    const std::vector<Case> cases = {
        // 257 - 1 = 2^8, B-smooth for B >= 256
        {257, 263, 512, 8},
        // 241 - 1 = 2^4 * 3 * 5
        {241, 547, 64, 8},
        // 101 - 1 = 2^2 * 5^2
        {101, 1009, 32, 8},
    };

    for (const auto& tc : cases) {
        const mpz_class p = tc.p;
        const mpz_class q = tc.q;
        const mpz_class N = p * q;
        const auto factors = primetools::PollardsPMinus1<mpz_class>(N, tc.B, tc.bases);
        ASSERT_TRUE(factors.has_value());
        EXPECT_TRUE(factors->first == p || factors->first == q);
        EXPECT_TRUE(factors->second == p || factors->second == q);
    }
}

TEST(Pollard, PollardsPMinus1_EdgeCases)
{
    EXPECT_FALSE(primetools::PollardsPMinus1<mpz_class>(mpz_class(0)).has_value());
    EXPECT_FALSE(primetools::PollardsPMinus1<mpz_class>(mpz_class(1)).has_value());

    // Prime input should never yield a non-trivial factor.
    EXPECT_FALSE(primetools::PollardsPMinus1<mpz_class>(mpz_class(101), /*B=*/128, /*Bases=*/16).has_value());

    // Even composite: base a=2 shares gcd 2 with N, so we should immediately return 2.
    const mpz_class N = mpz_class(2) * mpz_class(263);
    const auto factors = primetools::PollardsPMinus1<mpz_class>(N, /*B=*/128, /*Bases=*/1);
    ASSERT_TRUE(factors.has_value());
    EXPECT_TRUE(factors->first == 2 || factors->second == 2);
}

TEST(Pollard, PollardsPMinus1MT)
{
    const mpz_class p = 257;
    const mpz_class q = 263;
    const mpz_class N = p * q;

    const auto factors = primetools::PollardsPMinus1MT<mpz_class>(N, /*Threads=*/2, /*B=*/512, /*Bases=*/32);
    ASSERT_TRUE(factors.has_value());
    EXPECT_TRUE(factors->first == p || factors->first == q);
    EXPECT_TRUE(factors->second == p || factors->second == q);
}
