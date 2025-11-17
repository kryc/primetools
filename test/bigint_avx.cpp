#include <array>
#include <limits>
#include <vector>

#include <gmpxx.h>
#include <gtest/gtest.h>

#include "bigint_avx.hpp"

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

TEST(BigIntAVX, AddEquals)
{
    BigIntAVX<1024> big;

    // Initialize all lanes to zero
    {
        std::array<uint64_t, 16> vals = {0};
        std::span<const uint64_t> s(vals.data(), vals.size());
        big = s;
    }

    // Add 1 to all lanes
    big += 1;

    // Verify all lanes are now 1
    for (uint32_t lane = 0; lane < 16; ++lane) {
        EXPECT_EQ(big.lane_word(lane, 0), 1u);
        EXPECT_EQ(big.lane_word(lane, 1), 0u);
    }

    // Add another 1
    big += 1;
    // Verify all lanes are now 2
    for (uint32_t lane = 0; lane < 16; ++lane) {
        EXPECT_EQ(big.lane_word(lane, 0), 2u);
        EXPECT_EQ(big.lane_word(lane, 1), 0u);
    }

    big += 2;
    // Verify all lanes are now 4
    for (uint32_t lane = 0; lane < 16; ++lane) {
        EXPECT_EQ(big.lane_word(lane, 0), 4u);
        EXPECT_EQ(big.lane_word(lane, 1), 0u);
    }

    // Overflow the low limb in all lanes
    big += 0xFFFFFFFCull;
    // Verify all lanes are now 0 (with carry to next limb)
    for (uint32_t lane = 0; lane < 16; ++lane) {
        EXPECT_EQ(big.lane_word(lane, 0), 0u);
        EXPECT_EQ(big.lane_word(lane, 1), 1u);
    }
}

TEST(BigIntAVX, SubtractEquals)
{
    BigIntAVX<1024> big;

    big = 100;
    big -= 50;
    for (uint32_t lane = 0; lane < 16; ++lane) {
        EXPECT_EQ(big.lane_word(lane, 0), 50u);
        EXPECT_EQ(big.lane_word(lane, 1), 0u);
    }

    big = 0xFFFFFFFFull; // 2^32 - 1
    big += 1;
    EXPECT_EQ(big.lane_word(0, 0), 0u);
    EXPECT_EQ(big.lane_word(0, 1), 1u);
    big -= 1;
    EXPECT_EQ(big.lane_word(0, 0), 0xFFFFFFFFu);
    EXPECT_EQ(big.lane_word(0, 1), 0u);

    big = 0xFFFFFFFFFFFFFFFFull; // 2^64 - 1
    big += 1;
    EXPECT_EQ(big.lane_word(0, 0), 0u);
    EXPECT_EQ(big.lane_word(0, 1), 0u);
    EXPECT_EQ(big.lane_word(0, 2), 1u);
    big -= 1;
    EXPECT_EQ(big.lane_word(0, 0), 0xFFFFFFFFu);
    EXPECT_EQ(big.lane_word(0, 1), 0xFFFFFFFFu);
    EXPECT_EQ(big.lane_word(0, 2), 0u);

    BigIntAVX<1024> big2;
    big = 100;
    big2 = 50;
    big -= big2;
    for (uint32_t lane = 0; lane < 16; ++lane) {
        EXPECT_EQ(big.lane_word(lane, 0), 50u);
        EXPECT_EQ(big.lane_word(lane, 1), 0u);
    }

    big = 0xFFFFFFFFFFFFFFFFull; // 2^64 - 1
    big2 = 1;
    big -= big2;
    for (uint32_t lane = 0; lane < 16; ++lane) {
        EXPECT_EQ(big.lane_word(lane, 0), 0xFFFFFFFEu);
        EXPECT_EQ(big.lane_word(lane, 1), 0xFFFFFFFFu);
    }
}

TEST(BigIntAVX, LeftShift)
{
    BigIntAVX<1024> big;

    // Initialize all lanes to 1
    big = 1;

    // Shift left by 1
    big <<= 1;
    for (uint32_t lane = 0; lane < 16; ++lane) {
        EXPECT_EQ(big.lane_word(lane, 0), 2u);
        EXPECT_EQ(big.lane_word(lane, 1), 0u);
    }

    // Shift left by 31 more (total 32)
    big <<= 31;
    for (uint32_t lane = 0; lane < 16; ++lane) {
        EXPECT_EQ(big.lane_word(lane, 0), 0u);
        EXPECT_EQ(big.lane_word(lane, 1), 1u);
    }

    // Shift left by 32 more (total 64)
    big <<= 32;
    for (uint32_t lane = 0; lane < 16; ++lane) {
        EXPECT_EQ(big.lane_word(lane, 0), 0u);
        EXPECT_EQ(big.lane_word(lane, 1), 0u);
        EXPECT_EQ(big.lane_word(lane, 2), 1u);
    }
}

TEST(BigIntAVX, Divides)
{
    BigIntAVX<1024> big1, big2;
    // Use all 16 lanes with a mix of divisible and non-divisible pairs.

    // Case 1: all lanes divisible (each lane: divisor = 100, dividend = 100 * lane)
    {
        std::array<uint64_t, 16> dividends{};
        std::array<uint64_t, 16> divisors{};
        for (uint32_t lane = 0; lane < 16; ++lane) {
            dividends[lane] = static_cast<uint64_t>(100ull * (lane + 1));
            divisors[lane]  = 100ull;
        }
        std::span<const uint64_t> s_dividends(dividends.data(), dividends.size());
        std::span<const uint64_t> s_divisors(divisors.data(), divisors.size());
        big1 = s_divisors;   // divisor
        big2 = s_dividends;  // dividend
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
        big1 = s_divisors;   // divisor
        big2 = s_dividends;  // dividend
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
        big1 = s_divisors;   // divisor
        big2 = s_dividends;  // dividend
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