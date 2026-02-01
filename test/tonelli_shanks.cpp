#include <gtest/gtest.h>

#include <cstdint>

#include "tonelli_shanks.hpp"

TEST(TonelliShanks, FindsSquareRootModPrime)
{
    // 6^2 = 36 ≡ 10 (mod 13)
    const uint64_t n = 10;
    const uint64_t p = 13;

    const auto r = primetools::TonelliShanks<uint64_t>(n, p);
    ASSERT_TRUE(r.has_value());

    const uint64_t x = *r;
    EXPECT_EQ((x * x) % p, n % p);
}

TEST(TonelliShanks, ReturnsNulloptForNonResidue)
{
    // 2 is not a quadratic residue mod 11
    const uint64_t n = 2;
    const uint64_t p = 11;

    const auto r = primetools::TonelliShanks<uint64_t>(n, p);
    EXPECT_FALSE(r.has_value());
}

TEST(TonelliShanks, HandlesMod2)
{
    EXPECT_EQ(*primetools::TonelliShanks<uint64_t>(0, 2), 0u);
    EXPECT_EQ(*primetools::TonelliShanks<uint64_t>(1, 2), 1u);
    EXPECT_EQ(*primetools::TonelliShanks<uint64_t>(3, 2), 1u);
}
