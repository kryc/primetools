
#include <gtest/gtest.h>
#include <gmpxx.h>

#include "wheel.hpp"

namespace primetools {

TEST(Prime, WheelGeneration210)
{
    auto wheel = GetWheelGapsForModulus(210);
    EXPECT_FALSE(wheel.empty());
    // The wheel should start with 1 and end at 211
    mpz_class current = 1;
    for (auto gapword : wheel) {
        for (size_t i = 0; i < kGapsPerWord; ++i) {
            const uint64_t gap = gapword & kGapMask;
            gapword >>= kBitsPerWheelGap;
            current += gap;
        }
    }
    EXPECT_EQ(current, mpz_class(211));
}

TEST(Prime, WheelGeneration2310)
{
    auto wheel = GetWheelGapsForModulus(2310);
    EXPECT_FALSE(wheel.empty());
    // The wheel should start with 1 and end at 2311
    mpz_class current = 1;
    for (auto gapword : wheel) {
        for (size_t i = 0; i < kGapsPerWord; ++i) {
            const uint64_t gap = gapword & kGapMask;
            gapword >>= kBitsPerWheelGap;
            current += gap;
        }
    }
    EXPECT_EQ(current, mpz_class(2311));
}

TEST(Prime, WheelGeneration30030)
{
    auto wheel = GetWheelGapsForModulus(30030);
    EXPECT_FALSE(wheel.empty());
    // The wheel should start with 1 and end at 30031
    mpz_class current = 1;
    for (auto gapword : wheel) {
        for (size_t i = 0; i < kGapsPerWord; ++i) {
            const uint64_t gap = gapword & kGapMask;
            gapword >>= kBitsPerWheelGap;
            current += gap;
        }
    }
    EXPECT_EQ(current, mpz_class(30031));
}

TEST(Prime, WheelGeneration510510)
{
    auto wheel = GetWheelGapsForModulus(510510);
    EXPECT_FALSE(wheel.empty());
    // The wheel should start with 1 and end at 510511
    mpz_class current = 1;
    for (auto gapword : wheel) {
        for (size_t i = 0; i < kGapsPerWord; ++i) {
            const uint64_t gap = gapword & kGapMask;
            gapword >>= kBitsPerWheelGap;
            current += gap;
        }
    }
    EXPECT_EQ(current, mpz_class(510511));
}

TEST(Prime, WheelGeneration9699690)
{
    auto wheel = GetWheelGapsForModulus(9699690);
    EXPECT_FALSE(wheel.empty());
    // The wheel should start with 1 and end at 9699691
    mpz_class current = 1;
    for (auto gapword : wheel) {
        for (size_t i = 0; i < kGapsPerWord; ++i) {
            const uint64_t gap = gapword & kGapMask;
            gapword >>= kBitsPerWheelGap;
            current += gap;
        }
    }
    EXPECT_EQ(current, mpz_class(9699691));
}

// TEST(Prime, WheelGeneration223092870)
// {
//     auto wheel = GetWheelGapsForModulus(223092870);
//     EXPECT_FALSE(wheel.empty());
//     // The wheel should start with 1 and end at 223092871
//     mpz_class current = 1;
//     for (auto gapword : wheel) {
//         for (size_t i = 0; i < kGapsPerWord; ++i) {
//             const uint64_t gap = gapword & kGapMask;
//             gapword >>= kBitsPerWheelGap;
//             current += gap;
//         }
//     }
//     EXPECT_EQ(current, mpz_class(223092871));
// }

} // namespace primetools