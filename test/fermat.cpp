#include <array>
#include <limits>
#include <vector>

#include <gtest/gtest.h>

#include "fermat.hpp"

// 32-bit primes that are close together
std::vector<std::array<uint32_t, 2>> FermatTestCases = {
    { { 4294967291, 4294967279 } },
    { { 4294967231, 4294967211 } },
    { { 4294967189, 4294967171 } },
    { { 4294967161, 4294967143 } },
    { { 4294967113, 4294967097 } }
};

TEST(Fermat, StandardFermat)
{
    for (const auto& testCase : FermatTestCases) {
        const mpz_class p = testCase[0];
        const mpz_class q = testCase[1];
        const mpz_class N = p * q;

        const auto factors = primetools::FermatFactorisationAlgorithm1(N, 0, std::numeric_limits<size_t>::max());
        ASSERT_TRUE(factors.has_value());
        EXPECT_TRUE(factors->first == p || factors->first == q);
        EXPECT_TRUE(factors->second == p || factors->second == q);
    }
}

TEST(Fermat, FermatFactorisationAlgorithm2)
{
    for (const auto& testCase : FermatTestCases) {
        const mpz_class p = testCase[0];
        const mpz_class q = testCase[1];
        const mpz_class N = p * q;

        const auto factors = primetools::FermatFactorisationAlgorithm2(N, std::numeric_limits<size_t>::max());
        ASSERT_TRUE(factors.has_value());
        EXPECT_TRUE(factors->first == p || factors->first == q);
        EXPECT_TRUE(factors->second == p || factors->second == q);
    }
}

TEST(Fermat, ModifiedFermatFactorisation4)
{
    for (const auto& testCase : FermatTestCases) {
        const mpz_class p = testCase[0];
        const mpz_class q = testCase[1];
        const mpz_class N = p * q;

        const auto factors = primetools::ModifiedFermatFactorisation4(N, std::numeric_limits<size_t>::max());
        ASSERT_TRUE(factors.has_value());
        EXPECT_TRUE(factors->first == p || factors->first == q);
        EXPECT_TRUE(factors->second == p || factors->second == q);
    }
}

TEST(Fermat, FMod20Precomp)
{
    std::vector<std::array<uint32_t, 2>> testCases = {
        { { 151413473, 156367177 } }, // N mod 20 == 1
        { { 239246351, 174550133 } }, // N mod 20 == 3
        { { 240075973, 204127919 } }, // N mod 20 == 7
        { { 145980617, 149979937 } }, // N mod 20 == 9
        { { 183691033, 156659227 } }, // N mod 20 == 11
        { { 206155219, 224283167 } }, // N mod 20 == 13
        { { 211208359, 156869803 } }, // N mod 20 == 17
        { { 154892299, 163950821 } }  // N mod 20 == 19
    };

    for (const auto& testCase : testCases) {
        const mpz_class p = testCase[0];
        const mpz_class q = testCase[1];
        const mpz_class N = p * q;

        const auto factors = primetools::FMMod20Precomp(N, std::numeric_limits<size_t>::max());
        ASSERT_TRUE(factors.has_value());
        EXPECT_TRUE(factors->first == p || factors->first == q);
        EXPECT_TRUE(factors->second == p || factors->second == q);
    }
}