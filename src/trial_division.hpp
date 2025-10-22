#ifndef TRIAL_DIVISION_HPP
#define TRIAL_DIVISION_HPP

#include <array>
#include <iostream>
#include <optional>
#include <variant>
#include <vector>

#include <gmpxx.h>

#include "analyse.hpp"
#include "random.hpp"
#include "trial_division_data.hpp"
#include "util.hpp"

namespace primetools {

const std::vector<__uint128_t>
GenerateWheelGapsForModulus(
    const size_t Modulus,
    const size_t BitSize = 5,
    const PackingType Packing = PackingType::FastPack
);

template <typename T>
static inline
const std::tuple<const T, const T, const size_t>
GetUpperAndLowerBounds(
    const T& N,
    const size_t Modulus,
    const bool GuessSize,
    const size_t Bits,
    const T& RangeLower = 0,
    const T& RangeUpper = 0
)
{
    const size_t bits = Bits != 0 ? Bits : (GuessSize ? primetools::GuessSizeOfPrimeFactors(N, true) : primetools::bit_size(N));
    T lower_bound = GuessSize ? (T(1) << (bits - 1)) + RangeLower : RangeLower;
    T upper_bound = (RangeUpper == 0) 
        ? (GuessSize ? primetools::min(T(1) << bits, T(sqrt(N))) : T(sqrt(N))) 
        : (GuessSize ? lower_bound + RangeUpper : RangeUpper);

    // Step back lower_bound to be a multiple of Modulus
    if (lower_bound % Modulus != 0) {
        lower_bound -= (lower_bound % Modulus);
    }

    // Increment upper_bound to be a multiple of Modulus
    if (upper_bound % Modulus != 0) {
        upper_bound += (Modulus - (upper_bound % Modulus));
    }

    return std::make_tuple(lower_bound, upper_bound, bits);
}

template <typename T, const size_t Modulus, const size_t BitSize, typename AT, const size_t Count, const PackingType Packed>
static const std::optional<std::pair<T, T>>
TrialDivisionRange(
    const T& N,
    const T& StartValue,
    const T& EndValue,
    const std::span<const AT, Count> GapArray
)
{
    constexpr size_t WordBits = sizeof(AT) * 8;
    constexpr size_t GapsPerWord =  Packed == FastPack ? (WordBits - 1) / BitSize : WordBits / BitSize;
    constexpr uint64_t GapsMask = Packed == FastPack ? (1 << (BitSize + 1)) - 2 : (1 << BitSize) - 1;
    // const size_t ArraySize = constexpr Count == std::dynamic_extent ? GapArray.size() : Count;

    T candidate = StartValue;

    // Cycle forwards until the candidate is congruent to 1 modulo the specified Modulus
    // First, ensure the candidate is odd
    if ((candidate & 1) == 0) {
        candidate += 1;
    }

    // Then scan forwards until candidate % Modulus == 1
    while (primetools::modulo(candidate, Modulus) != 1 && candidate > 1) {
        // Check if we found a factor
        if (primetools::divides(N, candidate) && candidate > 1) {
            return std::make_pair(candidate, N / candidate);
        }
        candidate += 2;
    }

    // Special case for wheel as we can use an optimization to rotate
    // the gaps using a bit rotation. This also avoids a branch.
    if constexpr (Modulus == 30)
    {
        uint32_t gapword = WHEEL30GAPSUINT32;
        while (candidate <= EndValue)
        {
            if (primetools::divides(N, candidate) && candidate > 1) {
                return std::make_pair(candidate, N / candidate);
            }

            primetools::increment(candidate, gapword & 0xf);

            gapword = std::rotr(gapword, 4);
        }
    }
    else
    {
        while (candidate <= EndValue)
        {
            for (auto gapword : GapArray)
            {
                for (size_t index = 0; index < GapsPerWord; index++)
                {
                    if (primetools::divides(N, candidate) && candidate > 1) {
                        return std::make_pair(candidate, N / candidate);
                    }

                    if constexpr(Packed == DensePack) {
                        primetools::increment(candidate, (gapword & GapsMask) << 1);
                    } else {
                        primetools::increment(candidate, gapword & GapsMask);
                    }

                    gapword >>= BitSize;
                }
            }
        }
    }

    return std::nullopt;
}

template <typename T>
static const std::optional<std::pair<T, T>>
TrialDivisionRange(
    const T& N,
    const T& StartValue,
    const T& EndValue,
    const size_t Modulus = 510510
)
{
    switch(Modulus) {
        case 200560490130:
        {
            std::cout << "Generating wheel gaps for modulus " << Modulus << std::endl;
            auto gaps = GenerateWheelGapsForModulus(Modulus, 5, FastPack);
            std::cout << "Factorising using wheel" << Modulus << " factorization." << std::endl;
            return TrialDivisionRange<T, 200560490130, 5, __uint128_t, std::dynamic_extent, FastPack>(N, StartValue, EndValue, gaps);
        }
        case 6469693230:
        {
            std::cout << "Generating wheel gaps for modulus " << Modulus << std::endl;
            auto gaps = GenerateWheelGapsForModulus(Modulus, 5, FastPack);
            std::cout << "Factorising using wheel" << Modulus << " factorization." << std::endl;
            return TrialDivisionRange<T, 6469693230, 5, __uint128_t, std::dynamic_extent, FastPack>(N, StartValue, EndValue, gaps);
        }
        case 223092870:
            return TrialDivisionRange<T, 223092870, 5, __uint128_t, WHEEL223092870GAP_COUNT, FastPack>(N, StartValue, EndValue, WHEEL223092870GAPS);
        case 9699690:
            return TrialDivisionRange<T, 9699690, 5, __uint128_t, WHEEL9699690GAP_COUNT, FastPack>(N, StartValue, EndValue, WHEEL9699690GAPS);
        case 510510:
            return TrialDivisionRange<T, 510510, 4, uint64_t, WHEEL510510GAP_COUNT, DensePack>(N, StartValue, EndValue, WHEEL510510GAPS);
        case 30030:
            return TrialDivisionRange<T, 30030, 4, uint64_t, WHEEL30030GAP_COUNT, DensePack>(N, StartValue, EndValue, WHEEL30030GAPS);
        case 2310:
            return TrialDivisionRange<T, 2310, 4, uint64_t, WHEEL2310GAP_COUNT, Unpacked>(N, StartValue, EndValue, WHEEL2310GAPS);
        case 210:
            return TrialDivisionRange<T, 210, 4, uint64_t, WHEEL210GAP_COUNT, Unpacked>(N, StartValue, EndValue, WHEEL210GAPS);
        case 30:
            return TrialDivisionRange<T, 30, 4, uint64_t, WHEEL30GAP_COUNT, Unpacked>(N, StartValue, EndValue, WHEEL30GAPS);
        default:
            std::cerr << "Error: Unsupported modulus for trial division: " << Modulus << std::endl;
            return std::nullopt;
    }
}

template <typename T>
static const std::optional<std::pair<T, T>>
TrialDivision(
    const T& N,
    const size_t Modulus,
    const bool GuessSize = true,
    const size_t Bits = 0,
    const T RangeLower = 0,
    const T RangeUpper = 0
)
{
    if (N < 2) {
        return std::nullopt;
    }

    // Calculate the bounds and bits
    auto [lower_bound, upper_bound, bits] = GetUpperAndLowerBounds<T>(N, Modulus, GuessSize, Bits, RangeLower, RangeUpper);

    if (GuessSize) {
        std::cout << "Trying factorization of " << bits << "-bit primes using wheel" << Modulus << " factorization." << std::endl;
    }
    std::cout << "Searching primes in range [" << lower_bound << ", " << upper_bound <<
        "] using wheel" << Modulus << " factorization." << std::endl;

    return TrialDivisionRange<T>(N, lower_bound, upper_bound, Modulus);
}

template <typename T>
const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionBitflip(
    const T& N,
    const size_t MaxIterations,
    const bool UseGuessSize = true
)
{
    if (N < 2) {
        return std::nullopt;
    }

    // Calculate the size of N's constituent primes
    const size_t bits = primetools::GuessSizeOfPrimeFactors(N, true);
    
    // Get the square root of N as this is our upper bound for trial division
    const T upper_bound = T(sqrt(N));
    const size_t upper_bound_bits = primetools::bit_size(upper_bound);

    std::cout << "Trying random factorization of " << bits << "-bit primes " << std::endl;
    std::cout << "WARNING: This is highly inefficient and not recommended!" << std::endl;
    
    // We want to toggle bits between the highest and lowest bit
    const size_t bit_range = UseGuessSize ? bits - 2 : upper_bound_bits - 2; // Exclude highest and lowest bit

    // Make the first candidate the highest and lowest bit set
    T candidate = (T(1) << (bit_range + 1)) | 1;

    // Set up our PRNG
    primetools::MiniPRNG32 prng;

    for (size_t i = 0; i < MaxIterations; ++i) {
        // Flip a random bit in the candidate
        const size_t bit = (prng.Next() % bit_range) + 1;
        primetools::toggle_bit(candidate, bit);

        // Check if it divides N
        if (primetools::divides(N, candidate)) {
            return std::make_pair(candidate, N / candidate);
        }
    }

    return std::nullopt;
};

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionRandom(
    const mpz_class& N,
    const bool GuessSize = true,
    const size_t Bits = 0,
    const mpz_class& RangeLower = 0,
    const mpz_class& RangeUpper = 0,
    const uint64_t Seed = 0,
    const size_t Modulus = 510510,
    const size_t MaxIterations = std::numeric_limits<size_t>::max()
);

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionSimd(
    const mpz_class& N,
    const size_t MaxIterations
);

template <typename T>
inline std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionLinear(const T& N, const size_t Base, const size_t MaxIterations) {
    return TrialDivisionWheel510510<T>(N, MaxIterations);
}

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionMT(
    const mpz_class& N,
    const size_t Threads = 0,
    const size_t BlockSize = 0,
    const bool GuessSize = true,
    const size_t Bits = 0,
    const mpz_class& RangeLower = 0,
    const mpz_class& RangeUpper = 0,
    const size_t Modulus = 510510
);

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionRandomMT(
    const mpz_class& N,
    const size_t Threads = 0,
    const size_t BlockSize = 0,
    const bool GuessSize = true,
    const size_t Bits = 0,
    const mpz_class& RangeLower = 0,
    const mpz_class& RangeUpper = 0,
    const size_t Modulus = 510510
);

} // namespace primetools

#endif // TRIAL_DIVISION_HPP