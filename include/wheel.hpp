#ifndef WHEEL_HPP
#define WHEEL_HPP

#include <vector>

#include "prime.hpp"
#include "trial_division_data.hpp"

namespace primetools {

constexpr size_t kBitsPerWheelGap = 6;
constexpr uint64_t kGapMask = (1 << kBitsPerWheelGap) - 1;
constexpr size_t kMaxWheelGap = (1 << kBitsPerWheelGap) - 1;
constexpr size_t kGapsPerWord = (sizeof(uint64_t) * 8) / kBitsPerWheelGap;

constexpr uint32_t kWheel30 = 0x26424246;
constexpr uint32_t kWheel30BitsPerGap = 4;
constexpr uint32_t kWheel30Mask = (1 << kWheel30BitsPerGap) - 1;

const std::span<const uint64_t>
GetPrimesForWheelModulus(
    const size_t Modulus
);

const std::vector<uint64_t>
GenerateWheelGapsForModulus(
    const size_t Modulus
);

const std::span<const uint64_t>
GetWheelGapsForModulus(
    const size_t Modulus
);

} // namespace primetools

#endif // WHEEL_HPP