#include <gtest/gtest.h>

#include "primegenerator.hpp"

TEST(PrimeGenerator, PossibleSmallPrimes)
{
    primetools::PossiblePrimeGenerator<uint64_t> generator;
    std::vector<uint64_t> expected_primes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
        31, 37, 41, 43, 47, 53, 59, 61, 67, 71
    };
    for (const auto& prime : expected_primes) {
        EXPECT_EQ(generator.Next(), prime);
    }
}

TEST(PrimeGenerator, SmallPrimes)
{
    primetools::PrimeGenerator<uint64_t> generator;
    std::vector<uint64_t> expected_primes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
        31, 37, 41, 43, 47, 53, 59, 61, 67, 71
    };
    for (const auto& prime : expected_primes) {
        EXPECT_EQ(generator.Next(), prime);
    }
}

TEST(PrimeGenerator, LargePrimes210)
{
    primetools::PrimeGenerator<uint64_t> generator(210);
    // Advance to the 1000th prime
    for (size_t i = 0; i < 999; ++i) {
        generator.Next();
    }
    std::vector<uint64_t> expected_primes = {
        7919, 7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009
    };
    for (const auto& prime : expected_primes) {
        EXPECT_EQ(generator.Next(), prime);
    }
    // Advance to the millionth prime
    for (size_t i = 1009; i < 1'000'000; ++i) {
        generator.Next();
    }
    EXPECT_EQ(generator.Next(), 15'485'863);
}

TEST(PrimeGenerator, LargePrimes2310)
{
    primetools::PrimeGenerator<uint64_t> generator(2310);
    // Advance to the 1000th candidate
    for (size_t i = 0; i < 999; ++i) {
        generator.Next();
    }
    std::vector<uint64_t> expected_candidates = {
        7919, 7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009
    };
    for (const auto& candidate : expected_candidates) {
        EXPECT_EQ(generator.Next(), candidate);
    }
    // Advance to the millionth candidate
    for (size_t i = 1009; i < 1'000'000; ++i) {
        generator.Next();
    }
    EXPECT_EQ(generator.Next(), 15'485'863);
}

TEST(PrimeGenerator, LargePrimes30030)
{
    primetools::PrimeGenerator<uint64_t> generator(30030);
    // Advance to the 1000th candidate
    for (size_t i = 0; i < 999; ++i) {
        generator.Next();
    }
    std::vector<uint64_t> expected_candidates = {
        7919, 7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009
    };
    for (const auto& candidate : expected_candidates) {
        EXPECT_EQ(generator.Next(), candidate);
    }
    // Advance to the millionth candidate
    for (size_t i = 1009; i < 1'000'000; ++i) {
        generator.Next();
    }
    EXPECT_EQ(generator.Next(), 15'485'863);
}

TEST(PrimeGenerator, LargePrimes510510)
{
    primetools::PrimeGenerator<uint64_t> generator(510510);
    // Advance to the 1000th prime
    for (size_t i = 0; i < 999; ++i) {
        generator.Next();
    }
    std::vector<uint64_t> expected_primes = {
        7919, 7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009
    };
    for (const auto& prime : expected_primes) {
        EXPECT_EQ(generator.Next(), prime);
    }
    // Advance to the millionth prime
    for (size_t i = 1009; i < 1'000'000; ++i) {
        generator.Next();
    }
    EXPECT_EQ(generator.Next(), 15'485'863);
}