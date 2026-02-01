
#include <gtest/gtest.h>

#include "miller_rabin.hpp"

TEST(MillerRabin, SmallPrimes)
{
    std::vector<uint64_t> small_primes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
        31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
        73, 79, 83, 89, 97
    };

    for (const auto& prime : small_primes) {
        EXPECT_TRUE(primetools::IsPrime(prime)) << "Failed for prime: " << prime;
    }
}

TEST(MillerRabin, SmallComposites)
{
    std::vector<uint64_t> small_composites = {
        4, 6, 8, 9, 10, 12, 14, 15, 16, 18,
        20, 21, 22, 24, 25, 26, 27, 28, 30, 32,
        33, 34, 35, 36, 38, 39, 40, 42, 44, 45
    };

    for (const auto& composite : small_composites) {
        EXPECT_FALSE(primetools::IsPrime(composite)) << "Failed for composite: " << composite;
    }
}

TEST(MillerRabin, ProblemPrimes)
{
    std::vector<__uint128_t> problem_primes_u128 = {
        (static_cast<__uint128_t>(0xc) << 64) | 0x5a3c540f0c09b8f1  //227863093270133979377
    };

    for (const auto& prime : problem_primes_u128) {
        EXPECT_TRUE(primetools::IsPrime(prime)) << "Failed for prime: " << primetools::ToString(prime);
    }
}