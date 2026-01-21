#include <array>
#include <limits>
#include <vector>

#include <gtest/gtest.h>

#include "prime.hpp"
#include "trial_division.hpp"

// static const std::array<uint8_t, 8> ExpectedGaps30 = {
//     6, 4, 2, 4, 2, 4, 6, 2
// };
// static const std::array<uint8_t, 48> ExpectedGaps210 = {
//     10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4,
//     2, 6, 4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 6, 2, 4,
//     6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2
// };

// TEST(TrialDivision, GenerateWheelGapsForModulus30Unpacked)
// {
//     static constexpr size_t BitSize = 5;
//     __uint128_t expected = 0;
//     for (size_t i = 0; i < ExpectedGaps30.size(); i++) {
//         const size_t shift = i * BitSize;
//         expected |= (__uint128_t)(ExpectedGaps30[i]) << shift;
//     }
//     // Add the final gap to complete the wheel
//     auto gaps = primetools::GenerateWheelGapsForModulus(30, BitSize, primetools::PackingType::Unpacked);
//     ASSERT_EQ(gaps.size(), 1);
//     EXPECT_EQ(gaps.at(0), expected);
// }

// TEST(TrialDivision, GenerateWheelGapsForModulus30FastPack)
// {
//     static constexpr size_t BitSize = 5;
//     __uint128_t expected = 0;
//     for (size_t i = 0; i < ExpectedGaps30.size(); i++) {
//         const size_t shift = i * BitSize;
//         expected |= (__uint128_t)(ExpectedGaps30[i] >> 1) << shift;
//     }
//     expected <<= 1;
//     auto gaps = primetools::GenerateWheelGapsForModulus(30, BitSize, primetools::PackingType::FastPack);
//     ASSERT_EQ(gaps.size(), 1);
//     EXPECT_EQ(gaps.at(0), expected);
// }

// TEST(TrialDivision, GenerateWheelGapsForModulus210Unpacked)
// {
//     static constexpr size_t BitSize = 5;
//     static constexpr size_t WordBits = sizeof(__uint128_t) * 8;
//     static constexpr size_t GapsPerWord =  WordBits / BitSize;
//     std::vector<__uint128_t> expected;
//     __uint128_t current = 0;
//     size_t gap_count = 0;
//     for (size_t i = 0; i < ExpectedGaps210.size(); i++) {
//         const size_t shift = gap_count * BitSize;
//         current |= (__uint128_t)(ExpectedGaps210[i]) << shift;
//         gap_count++;
//         if (gap_count == GapsPerWord) {
//             expected.push_back(current);
//             current = 0;
//             gap_count = 0;
//         }
//     }
//     if (gap_count > 0) {
//         expected.push_back(current);
//     }
//     auto gaps = primetools::GenerateWheelGapsForModulus(210, BitSize, primetools::PackingType::Unpacked);
//     ASSERT_EQ(gaps.size(), expected.size());
//     for (size_t i = 0; i < gaps.size(); i++) {
//         EXPECT_EQ(gaps.at(i), expected.at(i));
//     }
// }

// TEST(TrialDivision, GenerateWheelGapsForModulus210FastPack)
// {
//     static constexpr size_t BitSize = 5;
//     static constexpr size_t WordBits = sizeof(__uint128_t) * 8;
//     static constexpr size_t GapsPerWord =  (WordBits - 1) / BitSize;
//     std::vector<__uint128_t> expected;
//     __uint128_t current = 0;
//     size_t gap_count = 0;
//     for (size_t i = 0; i < ExpectedGaps210.size(); i++) {
//         const size_t shift = gap_count * BitSize;
//         current |= (__uint128_t)(ExpectedGaps210[i] >> 1) << shift;
//         gap_count++;
//         if (gap_count == GapsPerWord) {
//             current <<= 1;
//             expected.push_back(current);
//             current = 0;
//             gap_count = 0;
//         }
//     }
//     if (gap_count > 0) {
//         current <<= 1;
//         expected.push_back(current);
//     }
//     auto gaps = primetools::GenerateWheelGapsForModulus(210, BitSize, primetools::PackingType::FastPack);
//     ASSERT_EQ(gaps.size(), expected.size());
//     for (size_t i = 0; i < gaps.size(); i++) {
//         EXPECT_EQ(gaps.at(i), expected.at(i));
//     }
// }

// TEST(TrialDivision, GenerateWheelGapsForModulus9699690)
// {
//     static constexpr size_t BitSize = 5;
//     auto gaps = primetools::GenerateWheelGapsForModulus(9699690, BitSize, primetools::PackingType::FastPack);
//     ASSERT_EQ(gaps.size(), primetools::WHEEL9699690GAPS.size());
//     for (size_t i =0; i < gaps.size(); i++) {
//         EXPECT_EQ(gaps[i], primetools::WHEEL9699690GAPS[i]);
//     }
// }

#ifdef PRIMETOOLS_ENABLE_LARGE_WHEEL_TEST
TEST(TrialDivision, GenerateWheelGapsForModulus223092870)
{
    static constexpr size_t BitSize = 5;
    auto gaps = primetools::GenerateWheelGapsForModulus(223092870, BitSize, primetools::PackingType::FastPack);
    ASSERT_EQ(gaps.size(), primetools::WHEEL223092870GAPS.size());
    for (size_t i =0; i < gaps.size(); i++) {
        EXPECT_EQ(gaps[i], primetools::WHEEL223092870GAPS[i]);
    }
}
#endif

TEST(TrialDivision, MeetInTheMiddleDoesNotSkipChunks)
{
    // Construct a semiprime where the smaller factor lies in a specific chunk.
    // This guards against gaps caused by incorrect meet-in-the-middle scheduling.
    static constexpr uint64_t p = 1009;
    static constexpr uint64_t q = 1'000'003;
    static constexpr uint64_t n = p * q;

    static constexpr size_t threads = 6;
    static constexpr size_t modulus = 30;
    static constexpr size_t block_size = 300; // multiple of modulus

    // Search a fixed range so the factor is inside a known chunk.
    static constexpr uint64_t range_lower = 0;
    static constexpr uint64_t range_upper = 2'100; // 7 chunks of 300

    auto result = primetools::TrialDivision<uint64_t>(
        n,
        threads,
        block_size,
        false,     // GuessSize
        0,         // Bits
        range_lower,
        range_upper,
        modulus,
        true,     // Status
        primetools::TrialDivisionStrategy::MeetInTheMiddle
    );

    ASSERT_TRUE(result.has_value());
    EXPECT_TRUE(result->HasFactor(p));
}