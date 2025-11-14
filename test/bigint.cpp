#include <array>
#include <limits>
#include <vector>

#include <gtest/gtest.h>

#include "bigint.hpp"

TEST(BigInt, Assign)
{
    BigInt<1024> bigInt;
    bigInt = 4294967295; // 2^32 - 1
    EXPECT_EQ(bigInt[0], 4294967295u);
    EXPECT_EQ(bigInt[1], 0u);

    mpz_class largeValue = (mpz_class(1) << 64) - 1; // 2^64 - 1
    bigInt = largeValue;
    EXPECT_EQ(bigInt[0], 4294967295u);
    EXPECT_EQ(bigInt[1], 4294967295u);
    EXPECT_EQ(bigInt[2], 0u);

    mpz_random(largeValue.get_mpz_t(), 1024);
    bigInt = largeValue;
    for (size_t i = 0; i < 32; ++i) {
        EXPECT_EQ(bigInt[i], static_cast<uint32_t>(largeValue.get_ui() & 0xFFFFFFFFu));
        largeValue >>= 32;
    }
}

TEST(BigInt, AddEquals)
{
    BigInt<1024> bigInt;
    bigInt += 1;
    EXPECT_EQ(bigInt[0], 1u);
    bigInt += 4294967295; // 2^32 - 1
    EXPECT_EQ(bigInt[0], 0u);
    EXPECT_EQ(bigInt[1], 1u);
}

TEST(BigInt, Divides)
{
    BigInt<1024> bigInt1, bigInt2;
    bigInt1 = 100;
    bigInt2 = 25;
    EXPECT_TRUE(bigInt1.divides(bigInt2));
    mpz_class a, b, mult;
    // Use proper bit-size random generation (previously mpz_random used limb count, producing oversized values)
    gmp_randstate_t rs; gmp_randinit_default(rs); gmp_randseed_ui(rs, 42u);
    mpz_urandomb(a.get_mpz_t(), rs, 512); // 512-bit random
    mpz_urandomb(b.get_mpz_t(), rs, 256); // 256-bit random
    mult = a * b;
    bigInt1 = mult;
    bigInt2 = a;
    EXPECT_TRUE(bigInt1.divides(bigInt2));
    bigInt2 += 1;
    EXPECT_FALSE(bigInt1.divides(bigInt2));
    gmp_randclear(rs);
}

TEST(BigInt, RestoringDivides)
{
    BigInt<1024> bigInt1, bigInt2;
    bigInt1 = 100;
    bigInt2 = 25;
    EXPECT_TRUE(bigInt1.restoring_divides(bigInt2));

    gmp_randstate_t rs; gmp_randinit_default(rs); gmp_randseed_ui(rs, 123u);
    mpz_class a, b, mult;
    mpz_urandomb(a.get_mpz_t(), rs, 512);
    mpz_urandomb(b.get_mpz_t(), rs, 256);
    mult = a * b;
    bigInt1 = mult;
    bigInt2 = a;
    EXPECT_TRUE(bigInt1.restoring_divides(bigInt2));
    bigInt2 += 1;
    EXPECT_FALSE(bigInt1.restoring_divides(bigInt2));
    gmp_randclear(rs);
}

#ifdef __AVX512F__
TEST(BigIntAVX, AssignU64)
{
    BigIntAVX<1024> big;

    // Fill all lanes with unique 64-bit values: lane L gets (L << 32) | (2^32-1).
    std::array<uint64_t, 16> vals;
    for (uint32_t lane = 0; lane < 16; ++lane) {
        vals[lane] = (static_cast<uint64_t>(lane) << 32) | 0xFFFFFFFFull;
    }
    std::span<const uint64_t> s(vals.data(), vals.size());
    big = s;

    // Check every lane directly via lane_word.
    for (uint32_t lane = 0; lane < 16; ++lane) {
        EXPECT_EQ(big.lane_word(lane, 0), 0xFFFFFFFFu);
        EXPECT_EQ(big.lane_word(lane, 1), lane);
    }
}

TEST(BigIntAVX, AssignMpz)
{
    BigIntAVX<1024> big;

    mpz_class v = (mpz_class(1) << 64) - 1; // 2^64 - 1
    std::span<const mpz_class> s(&v, 1);
    big = s;

    EXPECT_EQ(big.lane_word(0, 0), 0xFFFFFFFFu);
    EXPECT_EQ(big.lane_word(0, 1), 0xFFFFFFFFu);
    for (uint32_t lane = 1; lane < 16; ++lane) {
        EXPECT_EQ(big.lane_word(lane, 0), 0u);
        EXPECT_EQ(big.lane_word(lane, 1), 0u);
    }
}

TEST(BigIntAVX, Divides)
{
    BigIntAVX<1024> big1, big2;
    // Use all 16 lanes with a mix of divisible and non-divisible pairs.

    // Case 1: all lanes divisible (each lane: dividend = 100 * lane, divisor = 100)
    {
        std::array<uint64_t, 16> dividends{};
        std::array<uint64_t, 16> divisors{};
        for (uint32_t lane = 0; lane < 16; ++lane) {
            dividends[lane] = static_cast<uint64_t>(100ull * (lane + 1));
            divisors[lane]  = 100ull;
        }
        std::span<const uint64_t> s_dividends(dividends.data(), dividends.size());
        std::span<const uint64_t> s_divisors(divisors.data(), divisors.size());
        big1 = s_dividends;
        big2 = s_divisors;
        EXPECT_TRUE(big1.restoring_divides(big2));
    }

    // Case 2: none of the lanes divisible (add 1 to every dividend)
    {
        std::array<uint64_t, 16> dividends{};
        std::array<uint64_t, 16> divisors{};
        for (uint32_t lane = 0; lane < 16; ++lane) {
            dividends[lane] = static_cast<uint64_t>(100ull * (lane + 1) + 1ull);
            divisors[lane]  = 100ull;
        }
        std::span<const uint64_t> s_dividends(dividends.data(), dividends.size());
        std::span<const uint64_t> s_divisors(divisors.data(), divisors.size());
        big1 = s_dividends;
        big2 = s_divisors;
        EXPECT_FALSE(big1.restoring_divides(big2));
    }

    // Case 3: mixture of divisible and non-divisible lanes. Since restoring_divides
    // now returns true if any lane is divisible, expect true.
    {
        std::array<uint64_t, 16> dividends{};
        std::array<uint64_t, 16> divisors{};
        for (uint32_t lane = 0; lane < 16; ++lane) {
            bool divisible = (lane % 2 == 0);
            uint64_t base = static_cast<uint64_t>(10ull * (lane + 1));
            dividends[lane] = divisible ? base : (base + 1ull);
            divisors[lane]  = 10ull;
        }
        std::span<const uint64_t> s_dividends(dividends.data(), dividends.size());
        std::span<const uint64_t> s_divisors(divisors.data(), divisors.size());
        big1 = s_dividends;
        big2 = s_divisors;
        EXPECT_TRUE(big1.restoring_divides(big2));
    }
}

TEST(BigIntAVX, MaxBitLength)
{
    BigIntAVX<1024> big;

    // Test with all lanes zero
    {
        std::array<uint64_t, 16> vals = {0};
        std::span<const uint64_t> s(vals.data(), vals.size());
        big = s;
        EXPECT_EQ(big.max_bitlength(), 0u);
    }

    // Test with lane 0 set to 2^54
    {
        std::array<uint64_t, 16> vals = {0};
        vals[0] = uint64_t(1) << 54;
        std::span<const uint64_t> s(vals.data(), vals.size());
        big = s;
        EXPECT_EQ(big.max_bitlength(), 55u);
    }

    // Test with lane 15 set to max 1024-bit value
    {
        std::array<uint64_t, 16> vals = {0};
        vals[15] = std::numeric_limits<uint64_t>::max();
        std::span<const uint64_t> s(vals.data(), vals.size());
        big = s;
        EXPECT_EQ(big.max_bitlength(), 64);
    }

    // Test with lane 15 set to 1 << 1000
    {
        std::array<mpz_class, 16> vals = {0};
        vals[15] = mpz_class(1) << 1000;
        big = vals;
        EXPECT_EQ(big.max_bitlength(), 1001u);
    }
}


#endif