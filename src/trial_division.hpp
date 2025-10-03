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
using Range = std::pair<const T, const T>;

template <typename T>
static inline
const std::tuple<const T, const T, const size_t>
GetUpperAndLowerBounds(
    const T& N,
    const size_t Modulus,
    const bool GuessSize,
    const size_t Bits,
    const Range<T>& Range
)
{
    const size_t bits = Bits != 0 ? Bits : (GuessSize ? primetools::GuessSizeOfPrimeFactors(N, true) : primetools::bit_size(N));
    T lower_bound = GuessSize ? (T(1) << (bits - 1)) + Range.first : Range.first;
    T upper_bound = (Range.second == 0) 
        ? (GuessSize ? T(1) << bits : T(sqrt(N))) 
        : (GuessSize ? lower_bound + Range.second : Range.second);

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

template <typename T, const size_t Modulus, const size_t BitSize, const size_t Count, const std::array<const uint64_t, Count>& GapArray>
static const std::optional<std::pair<T, T>>
TrialDivisionRange(
    const T& N,
    const T& StartValue,
    const T& EndValue
)
{
    constexpr size_t GapsPerWord = 64 / BitSize;
    constexpr uint64_t GapsMask = (1 << BitSize) - 1;

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
            return TrialDivisionRange<T, 9699690, 6, 165888, WHEEL9699690GAPS>(N, StartValue, EndValue);
        case 510510:
            return TrialDivisionRange<T, 510510, 5, 7680, WHEEL510510GAPS>(N, StartValue, EndValue);
        case 30030:
            return TrialDivisionRange<T, 30030, 5, 480, WHEEL30030GAPS>(N, StartValue, EndValue);
        case 2310:
            return TrialDivisionRange<T, 2310, 4, 30, WHEEL2310GAPS>(N, StartValue, EndValue);
        case 210:
            return TrialDivisionRange<T, 210, 4, 3, WHEEL210GAPS>(N, StartValue, EndValue);
        case 30:
            return TrialDivisionRange<T, 30, 4, 1, WHEEL30GAPS>(N, StartValue, EndValue);
        default:
            std::cerr << "Error: Unsupported modulus for trial division: " << Modulus << std::endl;
            return std::nullopt;
    }
}

template <typename T, const size_t Modulus, const size_t BitSize, const size_t Count, const std::array<const uint64_t, Count>& GapArray>
static const std::optional<std::pair<T, T>>
TrialDivisionWheelInternal(
    const T& N,
    const bool GuessSize = true,
    const size_t Bits = 0,
    const Range<T>& Range = {T(0), T(0)}
)
{
    if (N < 2) {
        return std::nullopt;
    }

    // Calculate the bounds and bits
    auto [lower_bound, upper_bound, bits] = GetUpperAndLowerBounds<T>(N, Modulus, GuessSize, Bits, Range);

    if (GuessSize) {
        std::cout << "Trying factorization of " << bits << "-bit primes below " << upper_bound <<
            " using wheel" << Modulus << " factorization." << std::endl;
    }
    std::cout << "Trying factorization of primes in range [" << lower_bound << ", " << upper_bound <<
        "] using wheel" << Modulus << " factorization." << std::endl;

    return TrialDivisionRange<T, Modulus, BitSize, Count, GapArray>(N, lower_bound, upper_bound);
}

template <typename T>
static inline const std::optional<std::pair<T, T>>
TrialDivisionWheel9699690(const T& N, const bool GuessSize = true, const size_t Bits = 0, const T StartValue = 0, const T EndValue = 0) {
    Range<T> RangeValue = Range<T>{StartValue, EndValue};
    return TrialDivisionWheelInternal<T, 9699690, 6, 165888, WHEEL9699690GAPS>(N, GuessSize, Bits, RangeValue);
}

template <typename T>
static inline const std::optional<std::pair<T, T>>
TrialDivisionWheel510510(const T& N, const bool GuessSize = true, const size_t Bits = 0, const T StartValue = 0, const T EndValue = 0) {
    Range<T> RangeValue = Range<T>{StartValue, EndValue};
    return TrialDivisionWheelInternal<T, 510510, 5, 7680, WHEEL510510GAPS>(N, GuessSize, Bits, RangeValue);
}

template <typename T>
static inline const std::optional<std::pair<T, T>>
TrialDivisionWheel30030(const T& N, const bool GuessSize = true, const size_t Bits = 0, const T StartValue = 0, const T EndValue = 0) {
    Range<T> RangeValue = Range<T>{StartValue, EndValue};
    return TrialDivisionWheelInternal<T, 30030, 5, 480, WHEEL30030GAPS>(N, GuessSize, Bits, RangeValue);
}

template <typename T>
static inline const std::optional<std::pair<T, T>>
TrialDivisionWheel2310(const T& N, const bool GuessSize = true, const size_t Bits = 0, const T StartValue = 0, const T EndValue = 0) {
    Range<T> RangeValue = Range<T>{StartValue, EndValue};
    return TrialDivisionWheelInternal<T, 2310, 4, 30, WHEEL2310GAPS>(N, GuessSize, Bits, RangeValue);
}

template <typename T>
static inline const std::optional<std::pair<T, T>>
TrialDivisionWheel210(const T& N, const bool GuessSize = true, const size_t Bits = 0, const T StartValue = 0, const T EndValue = 0) {
    Range<T> RangeValue = Range<T>{StartValue, EndValue};
    return TrialDivisionWheelInternal<T, 210, 4, 3, WHEEL210GAPS>(N, GuessSize, Bits, RangeValue);
}

template <typename T>
static inline const std::optional<std::pair<T, T>>
TrialDivisionWheel30(const T& N, const bool GuessSize = true, const size_t Bits = 0, const T StartValue = 0, const T EndValue = 0) {
    Range<T> RangeValue = Range<T>{StartValue, EndValue};
    return TrialDivisionWheelInternal<T, 30, 4, 1, WHEEL30GAPS>(N, GuessSize, Bits, RangeValue);
}

template <typename T>
static const std::optional<std::pair<T, T>>
TrialDivision(
    const T& N,
    const size_t Modulus,
    const bool GuessSize = true,
    const size_t Bits = 0,
    const T StartValue = 0,
    const T EndValue = 0
)
{
    switch(Modulus) {
        case 9699690:
            return TrialDivisionWheel9699690<T>(N, GuessSize, Bits, StartValue, EndValue);
        case 510510:
            return TrialDivisionWheel510510<T>(N, GuessSize, Bits, StartValue, EndValue);
        case 30030:
            return TrialDivisionWheel30030<T>(N, GuessSize, Bits, StartValue, EndValue);
        case 2310:
            return TrialDivisionWheel2310<T>(N, GuessSize, Bits, StartValue, EndValue);
        case 210:
            return TrialDivisionWheel210<T>(N, GuessSize, Bits, StartValue, EndValue);
        case 30:
            return TrialDivisionWheel30<T>(N, GuessSize, Bits, StartValue, EndValue);
        default:
            std::cerr << "Error: Unsupported modulus for trial division: " << Modulus << std::endl;
            return std::nullopt;
    }
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

template <size_t BitSize = 5, size_t Modulus = 510510, size_t Count = 7680, const std::array<const uint64_t, Count>& GapArray = WHEEL510510GAPS>
const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionRandom(
    const mpz_class& N,
    const bool GuessSize = true,
    const size_t Bits = 0,
    const Range<mpz_class>& Range = {mpz_class(0), mpz_class(0)},
    const uint64_t Seed = 0,
    const size_t MaxIterations = std::numeric_limits<size_t>::max()
)
{
    constexpr size_t BlockSize = Modulus * 512;

    if (N < 2) {
        return std::nullopt;
    }

    // Calculate the bounds and bits
    auto [lower_bound, upper_bound, bits] = GetUpperAndLowerBounds<mpz_class>(N, Modulus, GuessSize, Bits, Range);

    // Calculate the difference between the bounds
    mpz_class diff = upper_bound - lower_bound;

    // Calculate the number of blocks
    mpz_class blocks = diff / BlockSize;

    std::cout << "Trying random factorization of primes of size " << bits << " bits. " << blocks << " blocks." << std::endl;
    std::cout << "WARNING: This is highly inefficient and not recommended!" << std::endl;

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, Seed);

    mpz_class block;
    mpz_urandomm(block.get_mpz_t(), state, blocks.get_mpz_t());

    for (size_t i = 0; i < MaxIterations; i++) {
        mpz_class block_start = lower_bound + (block * BlockSize);
        mpz_class block_end = block_start + BlockSize;

        auto result = TrialDivisionRange<mpz_class, Modulus, BitSize, Count, GapArray>(N, block_start, block_end);
        if (result) {
            return result;
        }
        
        // Move to a new random block
        mpz_urandomm(block.get_mpz_t(), state, blocks.get_mpz_t());
    }

    return std::nullopt;
}

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
    const bool GuessSize = true,
    const size_t Bits = 0,
    const Range<mpz_class>& Range = {mpz_class(0), mpz_class(0)},
    const size_t Modulus = 510510
);

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionRandomMT(
    const mpz_class& N,
    const size_t Threads = 0,
    const bool GuessSize = true,
    const size_t Bits = 0,
    const Range<mpz_class>& Range = {mpz_class(0), mpz_class(0)},
    const size_t Modulus = 510510
);

} // namespace primetools

#endif // TRIAL_DIVISION_HPP