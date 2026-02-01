#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <vector>

#include <gmpxx.h>

#include "prime.hpp"
#include "wheel.hpp"

namespace primetools {

namespace {
    // Mutex for getting generated modulus wheels
    std::mutex gWheelMutex;
    std::map<size_t, std::vector<uint64_t>> gWheelCache;
} // anonymous namespace

const std::span<const uint64_t>
GetWheelGapsForModulusLocked(
    const size_t Modulus
);

const std::span<const uint64_t>
GetPrimesForWheelModulus(
    const size_t Modulus
)
{
    std::span<const uint64_t> primes = GetSmallPrimes();
    switch(Modulus)
    {
        case 30:
            return primes.subspan(0, 3); // 2 - 5
        case 210:
            return primes.subspan(0, 4); // 2 - 7
        case 2310:
            return primes.subspan(0, 5); // 2 - 11
        case 30030:
            return primes.subspan(0, 6); // 2 - 13
        case 510510:
            return primes.subspan(0, 7); // 2 - 17
        case 9699690:
            return primes.subspan(0, 8); // 2 - 19
        case 223092870:
            return primes.subspan(0, 9); // 2 - 23
        case 6469693230:
            return primes.subspan(0, 10); // 2 - 29
        case 200560490130:
            return primes.subspan(0, 11); // 2 - 31
        case 7420738134810:
            return primes.subspan(0, 12); // 2 - 37
        case 304250263527210:
            return primes.subspan(0, 13); // 2 - 41
        default:
            throw std::invalid_argument("Unsupported wheel Modulus");
    }
}

const uint64_t
GetWheelStartValueForModulus(
    const size_t Modulus
)
{
    switch(Modulus)
    {
        case 30:
            return 7;
        case 210:
            return 11;
        case 2310:
            return 13;
        case 30030:
            return 17;
        case 510510:
            return 19;
        case 9699690:
            return 23;
        case 223092870:
            return 29;
        case 6469693230:
            return 31;
        case 200560490130:
            return 37;
        case 7420738134810:
            return 41;
        case 304250263527210:
            return 43;
        default:
            throw std::invalid_argument("Unsupported wheel Modulus for start value");
    }
}

inline static const uint64_t
GetPreviousModulus(
    const size_t Modulus
)
{
    switch(Modulus)
    {
        case 210:
            return 30;
        case 2310:
            return 210;
        case 30030:
            return 2310;
        case 510510:
            return 30030;
        case 9699690:
            return 510510;
        case 223092870:
            return 9699690;
        case 6469693230:
            return 223092870;
        case 200560490130:
            return 6469693230;
        case 7420738134810:
            return 200560490130;
        case 304250263527210:
            return 7420738134810;
        default:
            throw std::invalid_argument("Unsupported wheel Modulus for previous modulus");
    }
}

const std::vector<uint64_t>
GenerateWheelGapsForModulusSmallWheel(
    const size_t Modulus
)
{
    // Generate the gaps for the wheel
    std::vector<uint64_t> gaps;
    uint64_t next = 0;
    size_t count = 0;
    uint64_t last_residue = 1;
    uint32_t wheel30_gaps = std::rotr(kWheel30, kWheel30BitsPerGap);
    for (size_t residue = 7;
        residue < Modulus;
        residue += (wheel30_gaps & kWheel30Mask),
        wheel30_gaps = std::rotr(wheel30_gaps, kWheel30BitsPerGap))
    {
        if (std::gcd(residue, Modulus) != 1) {
            continue;
        }

        const uint64_t gap = residue - last_residue;
        if (gap > kMaxWheelGap) {
            throw std::runtime_error("Wheel gap exceeds maximum representable size.");
        }
        
        const size_t shift = count++ * kBitsPerWheelGap;
        next |= (gap << shift);
        if (count == kGapsPerWord) {
            gaps.push_back(next);
            next = 0;
            count = 0;
        }
        last_residue = residue;
    }

    // Add final gap to complete the cycle
    const uint64_t gap = Modulus - last_residue + 1;
    if (count > 0) {
        const size_t shift = count++ * kBitsPerWheelGap;
        next |= (gap << shift);
    } else {
        next = gap;
    }

    gaps.push_back(next);
    return gaps;
}

const std::vector<uint64_t>
GenerateWheelGapsForModulusLargeWheel(
    const size_t Modulus
)
{
    // Generate the gaps for the wheel
    std::vector<uint64_t> gaps;
    uint64_t next = 0;
    size_t count = 0;
    uint64_t last_residue = 1;

    // Generate the smaller wheel to use in constructing the gaps
    auto smaller_wheel = GetWheelGapsForModulusLocked(GetPreviousModulus(Modulus));
    uint64_t residue = 1;

    while (residue < Modulus) {
        for (uint64_t gapword : smaller_wheel) {
            for (size_t i = 0; i < kGapsPerWord; ++i) {
                const uint64_t nextStep = gapword & kGapMask;
                gapword >>= kBitsPerWheelGap;
                residue += nextStep;

                // std::cout << "Residue: " << residue << std::endl;

                // It can't exceed the modulus as the smaller wheels
                // are always a multiple of the larger wheel
                if (residue >= Modulus) {
                    break;
                }

                if (std::gcd(residue, Modulus) != 1) {
                    continue;
                }

                const uint64_t gap = residue - last_residue;
                if (gap > kMaxWheelGap) {
                    throw std::runtime_error("Wheel gap (" + std::to_string(gap) + ") exceeds maximum representable size.");
                }
                
                const size_t shift = count++ * kBitsPerWheelGap;
                next |= (gap << shift);
                if (count == kGapsPerWord) {
                    gaps.push_back(next);
                    next = 0;
                    count = 0;
                }
                last_residue = residue;
            }
        }
    }

    // Add final gap to complete the cycle
    const uint64_t gap = Modulus - last_residue + 1;
    if (count > 0) {
        const size_t shift = count++ * kBitsPerWheelGap;
        next |= (gap << shift);
    } else {
        next = gap;
    }

    gaps.push_back(next);
    return gaps;
}

const std::vector<uint64_t>
GenerateWheelGapsForModulus(
    const size_t Modulus
)
{
    if (Modulus <= 510510) {
        return GenerateWheelGapsForModulusSmallWheel(Modulus);
    } else {
        return GenerateWheelGapsForModulusLargeWheel(Modulus);
    }
}

const std::span<const uint64_t>
GetWheelGapsForModulusLocked(
    const size_t Modulus
)
{
    auto it = gWheelCache.find(Modulus);
    if (it != gWheelCache.end()) {
        return it->second;
    }
    // std::cout << "Generating wheel gaps for modulus " << Modulus << std::endl;
    auto gaps = GenerateWheelGapsForModulus(Modulus);
    // std::cout << "Generated " << gaps.size() << " gap words for wheel modulus " << Modulus << std::endl;
    gWheelCache[Modulus] = std::move(gaps);
    return gWheelCache[Modulus];
}

const std::span<const uint64_t>
GetWheelGapsForModulus(
    const size_t Modulus
)
{
    std::lock_guard<std::mutex> lock(gWheelMutex);
    return GetWheelGapsForModulusLocked(Modulus);
}

} // namespace primetools