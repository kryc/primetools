#include <array>
#include <limits>
#include <vector>

#include <gtest/gtest.h>

#include "prime.hpp"
#include "trial_division.hpp"
#include "trial_division_data.hpp"

TEST(TrialDivision, GenerateWheelGapsForModulus30)
{
    std::vector<uint8_t> expected_gaps = {
        6, 4, 2, 4, 2, 4, 6, 2
    };
    __uint128_t expected_gaps_128 = 0;
    for (size_t i = 0; i < expected_gaps.size(); i++) {
        const size_t shift = i * 5;
        expected_gaps_128 |= (__uint128_t)(expected_gaps[i]) << shift;
    }
    // Add the final gap to complete the wheel

    std::vector<__uint128_t> gaps = primetools::GenerateWheelGapsForModulus(30, 5, primetools::PackingType::FastPack);
    ASSERT_EQ(gaps.size(), 1);
    EXPECT_EQ(gaps.at(0), expected_gaps_128);
}

TEST(TrialDivision, GenerateWheelGapsForModulus210)
{
    std::vector<uint8_t> expected_gaps = {
        10, 2, 4, 2, 4, 6, 2, 6, 4, 2, 4, 6, 6, 2, 6, 4,
        2, 6, 4, 6, 8, 4, 2, 4, 2, 4, 8, 6, 4, 6, 2, 4,
        6, 2, 6, 6, 4, 2, 4, 6, 2, 6, 4, 2, 4, 2, 10, 2
    };
    std::vector<__uint128_t> expected_gaps_vec;
    __uint128_t current_gap_word = 0;
    size_t gap_count = 0;
    for (size_t i = 0; i < expected_gaps.size(); i++) {
        const size_t shift = gap_count * 5;
        current_gap_word |= (__uint128_t)(expected_gaps[i]) << shift;
        gap_count++;
        if (gap_count == 127 / 5) {
            expected_gaps_vec.push_back(current_gap_word);
            current_gap_word = 0;
            gap_count = 0;
        }
    }
    if (gap_count > 0) {
        expected_gaps_vec.push_back(current_gap_word);
    }
    std::vector<__uint128_t> gaps = primetools::GenerateWheelGapsForModulus(210, 5, primetools::PackingType::FastPack);
    ASSERT_EQ(gaps.size(), expected_gaps_vec.size());
    for (size_t i = 0; i < gaps.size(); i++) {
        EXPECT_EQ(gaps.at(i), expected_gaps_vec.at(i));
    }
}

TEST(TrialDivision, GenerateWheelGapsForModulus9699690)
{
    auto gaps = primetools::GenerateWheelGapsForModulus(9699690, 5, primetools::PackingType::FastPack);
    ASSERT_EQ(gaps.size(), primetools::WHEEL9699690GAPS.size());
    for (size_t i =0; i < gaps.size(); i++) {
        EXPECT_EQ(gaps[i], primetools::WHEEL9699690GAPS[i]);
    }
}

TEST(TrialDivision, GenerateWheelGapsForModulus223092870)
{
    auto gaps = primetools::GenerateWheelGapsForModulus(223092870, 5, primetools::PackingType::FastPack);
    ASSERT_EQ(gaps.size(), primetools::WHEEL223092870GAPS.size());
    for (size_t i =0; i < gaps.size(); i++) {
        EXPECT_EQ(gaps[i], primetools::WHEEL223092870GAPS[i]);
    }
}