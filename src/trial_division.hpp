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
using RangeOrGuess = std::variant<const Range<T>, bool>;

template <typename T>
const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel30(
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
    const T lower_bound = UseGuessSize ? T(1) << (bits - 1) : T(1);
    const T upper_bound = UseGuessSize ? (T(1) << bits) : T(sqrt(N));

    T candidate = lower_bound;

    // Cycle backwards until the candidate is congruent to 1 modulo 2310
    while (candidate % 30 != 1 && candidate > 1) {
        candidate -= 1;
    }

    std::cout << "Trying factorization of " << bits << "-bit primes below " << upper_bound <<
        " using wheel30 factorization." << std::endl << "Starting candidate: " << candidate << std::endl;

    uint32_t gapword = WHEEL30GAPS;

    for (size_t i = 0; i < MaxIterations && candidate <= upper_bound; ++i)
    {
        if (primetools::divides(N, candidate) && candidate > 1) {
            return std::make_pair(candidate, N / candidate);
        }

        // candidate += gap & 0xf;
        primetools::increment(candidate, gapword & 0xf);

        // Rotate the gaps
        gapword = std::rotr(gapword, 4); // rotate by 4 bits
    }

    return std::nullopt;
}

template <typename T, const size_t Modulus, const size_t BitSize, const size_t Count, const std::array<const uint64_t, Count>& GapArray>
const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel(
    const T& N,
    const RangeOrGuess<T>& RangeOrGuessSize = true
)
{
    constexpr size_t GapsPerWord = 64 / BitSize;
    constexpr uint64_t GapsMask = (1 << BitSize) - 1;
    T lower_bound = T(1);
    T upper_bound = T(sqrt(N));
    size_t bits = primetools::bit_size(N);

    if (N < 2) {
        return std::nullopt;
    }

    if (std::holds_alternative<bool>(RangeOrGuessSize) && std::get<bool>(RangeOrGuessSize)) {
        // Estimate the size of N's constituent primes
        bits = primetools::GuessSizeOfPrimeFactors(N, true);
        lower_bound = T(1) << (bits - 1);
        upper_bound = T(1) << bits;
    }
    else if (std::holds_alternative<const Range<T>>(RangeOrGuessSize)){
        lower_bound = std::get<const Range<T>>(RangeOrGuessSize).first;
        upper_bound = std::get<const Range<T>>(RangeOrGuessSize).second;
        if (lower_bound < 1) {
            lower_bound = T(1);
        }
        if (upper_bound < lower_bound) {
            upper_bound = T(sqrt(N));
        }
    }

    T candidate = lower_bound;

    // Cycle forwards until the candidate is congruent to 1 modulo the specified Modulus
    // First, ensure the candidate is odd
    if ((candidate & 1) == 0) {
        candidate += 1;
    }
    // Then scan forwards until candidate % Modulus == 1
    while (candidate % Modulus != 1 && candidate > 1) {
        // Check if we found a factor
        if (primetools::divides(N, candidate) && candidate > 1) {
            return std::make_pair(candidate, N / candidate);
        }
        candidate += 2;
    }

    std::cout << "Trying factorization of " << bits << "-bit primes below " << upper_bound <<
        " using wheel" << Modulus << " factorization." << std::endl;
    // std::cout << "Starting candidate: " << candidate << std::endl;

    while (candidate <= upper_bound)
    {
        for (auto gapword : GapArray)
        {
            for (size_t index = 0; index < GapsPerWord; index++)
            {
                if (primetools::divides(N, candidate) && candidate > 1) {
                    return std::make_pair(candidate, N / candidate);
                }

                // candidate += gap & 0x1f;
                primetools::increment(candidate, gapword & GapsMask);

                // Rotate the gaps
                gapword >>= BitSize; // shift by 5 bits
            }
        }
    }

    return std::nullopt;
}

template <typename T>
static inline const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel510510(const T& N, const bool GuessSize = true) {
    RangeOrGuess<T> RangeOrGuessSize = GuessSize;
    return TrialDivisionWheel<T, 510510, 5, 7680, WHEEL510510GAPS>(N, RangeOrGuessSize);
}
template <typename T>
static inline const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel510510(const T& N, const T StartValue, const T EndValue) {
    RangeOrGuess<T> RangeOrGuessSize = Range<T>{StartValue, EndValue};
    return TrialDivisionWheel<T, 510510, 5, 7680, WHEEL510510GAPS>(N, RangeOrGuessSize);
}

template <typename T>
static inline const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel30030(const T& N, const bool GuessSize = true) {
    RangeOrGuess<T> RangeOrGuessSize = GuessSize;
    return TrialDivisionWheel<T, 30030, 5, 480, WHEEL30030GAPS>(N, RangeOrGuessSize);
}
template <typename T>
static inline const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel30030(const T& N, const T StartValue, const T EndValue) {
    RangeOrGuess<T> RangeOrGuessSize = Range<T>{StartValue, EndValue};
    return TrialDivisionWheel<T, 30030, 5, 480, WHEEL30030GAPS>(N, RangeOrGuessSize);
}

template <typename T>
static inline const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel2310(const T& N, const bool GuessSize = true) {
    RangeOrGuess<T> RangeOrGuessSize = GuessSize;
    return TrialDivisionWheel<T, 2310, 4, 30, WHEEL2310GAPS>(N, RangeOrGuessSize);
}
template <typename T>
static inline const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel2310(const T& N, const T StartValue, const T EndValue) {
    RangeOrGuess<T> RangeOrGuessSize = Range<T>{StartValue, EndValue};
    return TrialDivisionWheel<T, 2310, 4, 30, WHEEL2310GAPS>(N, RangeOrGuessSize);
}

template <typename T>
static inline const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel210(const T& N, const bool GuessSize = true) {
    RangeOrGuess<T> RangeOrGuessSize = GuessSize;
    return TrialDivisionWheel<T, 210, 4, 3, WHEEL210GAPS>(N, RangeOrGuessSize);
}
template <typename T>
static inline const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel210(const T& N, const T StartValue, const T EndValue) {
    RangeOrGuess<T> RangeOrGuessSize = Range<T>{StartValue, EndValue};
    return TrialDivisionWheel<T, 210, 4, 3, WHEEL210GAPS>(N, RangeOrGuessSize);
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
TrialDivisionSimd(
    const mpz_class& N,
    const size_t MaxIterations
);

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionRandom(
    const mpz_class& N,
    const size_t Base,
    const size_t MaxIterations
);

template <typename T>
inline std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionLinear(const T& N, const size_t Base, const size_t MaxIterations) {
    return TrialDivisionWheel510510<T>(N, MaxIterations);
}

} // namespace primetools

#endif // TRIAL_DIVISION_HPP