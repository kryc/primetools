#include <array>
#include <limits>
#include <vector>

#include <gtest/gtest.h>

#include "prime.hpp"

TEST(Prime, GetNthPrime)
{
    std::vector<std::pair<uint64_t, uint64_t>> testCases = {
        { 1, 2 },
        { 2, 3 },
        { 3, 5 },
        { 4, 7 },
        { 5, 11 },
        { 6, 13 },
        { 7, 17 },
        { 10, 29 },
        { 25, 97 },
        { 100, 541 },
        { 1000, 7919 },
        { 10000, 104729 },
        { 100000, 1299709 }
    };

    for (const auto& testCase : testCases) {
        const uint64_t N = testCase.first;
        const uint64_t expectedPrime = testCase.second;

        const uint64_t prime = primetools::GetNthPrime(N);
        EXPECT_EQ(prime, expectedPrime);
    }
}

TEST(Prime, GetPrimesInRange)
{
    struct TestCase {
        uint64_t lower;
        uint64_t upper;
        std::vector<uint64_t> expectedPrimes;
    };

    std::vector<TestCase> testCases = {
        { 1, 1, { } },
        { 1, 2, { 2 } },
        { 1, 3, { 2, 3 } },
        { 1, 10, { 2, 3, 5, 7 } },
        { 10, 30, { 11, 13, 17, 19, 23, 29 } },
        { 50, 100, { 53, 59, 61, 67, 71, 73, 79, 83, 89, 97 } },
        { 100, 150, { 101, 103, 107, 109, 113, 127, 131, 137, 139, 149 } }
    };

    for (const auto& testCase : testCases) {
        const auto primes = primetools::GetPrimesInRange(testCase.lower, testCase.upper);
        EXPECT_EQ(primes, testCase.expectedPrimes);
    }
}

TEST(Prime, GetPrimesTo)
{
    struct TestCase {
        uint64_t upper;
        std::vector<uint64_t> expectedPrimes;
    };

    std::vector<TestCase> testCases = {
        { 3, { 2, 3 } },
        { 10, { 2, 3, 5, 7 } },
        { 30, { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 } },
        { 100, { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                 31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
                 73, 79, 83, 89, 97 } }
    };

    auto to1 = primetools::GetPrimesTo(1);
    EXPECT_EQ(to1.size(), 0);
    auto to2 = primetools::GetPrimesTo(2);
    EXPECT_EQ(to2.size(), 1);
    EXPECT_EQ(to2[0], 2);

    for (const auto& testCase : testCases) {
        const auto primes = primetools::GetPrimesTo(testCase.upper);
        EXPECT_EQ(primes.size(), testCase.expectedPrimes.size());
        for (size_t i = 0; i < primes.size(); i++) {
            EXPECT_EQ(primes[i], testCase.expectedPrimes[i]);
        }
    }
}

TEST(Prime, GetFirstNPrimes)
{
    struct TestCase {
        uint64_t N;
        std::vector<uint64_t> expectedPrimes;
    };

    std::vector<TestCase> testCases = {
        { 1, { 2 } },
        { 2, { 2, 3 } },
        { 3, { 2, 3, 5 } },
        { 5, { 2, 3, 5, 7, 11 } },
        { 10, { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29 } },
        { 20, { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
                31, 37, 41, 43, 47, 53, 59, 61, 67, 71 } }
    };

    for (const auto& testCase : testCases) {
        const auto primes = primetools::GetFirstNPrimes(testCase.N);
        EXPECT_EQ(primes, testCase.expectedPrimes);
    }
}