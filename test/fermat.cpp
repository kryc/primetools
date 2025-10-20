#include <array>
#include <limits>
#include <vector>

#include <gtest/gtest.h>

#include "fermat.hpp"

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
        EXPECT_EQ((factors->first == p && factors->second == q) ||
                  (factors->first == q && factors->second == p), true);
    }
}