#ifndef TRIAL_DIVISION_HPP
#define TRIAL_DIVISION_HPP

#include <array>
#include <iostream>
#include <optional>
#include <variant>

#include <gmpxx.h>

#include "analyse.hpp"
#include "random.hpp"
#include "trial_division_data.hpp"
#include "util.hpp"

namespace primetools {

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
        ? (GuessSize ? T(1) << bits : T(sqrt(N))) 
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

template <typename T, const size_t Modulus, const size_t BitSize, typename AT, const size_t Count, const bool Packed, const std::array<const AT, Count>& GapArray>
static const std::optional<std::pair<T, T>>
TrialDivisionRange(
    const T& N,
    const T& StartValue,
    const T& EndValue
)
{
    constexpr size_t GapsPerWord =  Packed ? 63 / BitSize : 64 / BitSize;
    constexpr uint64_t GapsMask = Packed ? (1 << (BitSize + 1)) - 2 : (1 << BitSize) - 1;

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

                    primetools::increment(candidate, gapword & GapsMask);

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
        case 9699690:
            return TrialDivisionRange<T, 9699690, 5, __uint128_t, 66356, true, WHEEL9699690GAPS>(N, StartValue, EndValue);
        case 510510:
            return TrialDivisionRange<T, 510510, 4, uint64_t, 6144, true, WHEEL510510GAPS>(N, StartValue, EndValue);
        case 30030:
            return TrialDivisionRange<T, 30030, 4, uint64_t, 384, true, WHEEL30030GAPS>(N, StartValue, EndValue);
        case 2310:
            return TrialDivisionRange<T, 2310, 4, uint64_t, 30, false, WHEEL2310GAPS>(N, StartValue, EndValue);
        case 210:
            return TrialDivisionRange<T, 210, 4, uint64_t, 3, false, WHEEL210GAPS>(N, StartValue, EndValue);
        case 30:
            return TrialDivisionRange<T, 30, 4, uint64_t, 1, false, WHEEL30GAPS>(N, StartValue, EndValue);
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