#ifndef TRIAL_DIVISION_HPP
#define TRIAL_DIVISION_HPP

#include <iostream>
#include <optional>

#include <gmpxx.h>

#include "analyse.hpp"
#include "random.hpp"
#include "util.hpp"

namespace primetools {

static const uint32_t WHEEL30GAPS = 0x26424246;
// 210‑wheel gaps packed 16×4‑bit per 64‑bit word
static const std::array<const uint64_t, 3> WHEEL210GAPS = {
    0x462664246264242aULL, // gaps[0..15]
    0x4264684242486462ULL, // gaps[16..31]
    0x2a24246264246626ULL  // gaps[32..47]
};
// 2310‑wheel gaps packed 16×4‑bit per 64‑bit word
static const std::array<const uint64_t, 30> WHEEL2310GAPS = {
    0x246266424626424cULL,
    0x62a264e424248646ULL,
    0x4242ac242a264246ULL,
    0x4864624626664626ULL,
    0x6424662642a68642ULL,
    0x24864242a2a62462ULL,
    0x684242c6462462c4ULL,
    0x42462642a6264264ULL,
    0x626646626642a2a2ULL,
    0x4264684624864624ULL,
    0x2a2a242462a24686ULL,
    0x6484626642468424ULL,
    0x4662666468424248ULL,
    0x6462a2a242462642ULL,
    0x4248a62486642462ULL,
    0x626424668426a842ULL,
    0x4246264242a2a264ULL,
    0x6842424864666266ULL,
    0x2424864246626484ULL,
    0x468642a264242a2aULL,
    0x6426468426486462ULL,
    0x42a2a24662664662ULL,
    0x64624626a2462642ULL,
    0x24c2642646c24248ULL,
    0x626426a2a2424684ULL,
    0x424686a246266424ULL,
    0x4626466626426468ULL,
    0x6642462a242ca242ULL,
    0x2646842424e462a2ULL,
    0x2c42462642466264ULL,
};

template <typename T>
const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel30(
    const T& N,
    const size_t MaxIterations
)
{
    if (N < 2) {
        return std::nullopt;
    }
    
    // Get the square root of N as this is our upper bound for trial division
    const T upper = sqrt(N);

    std::cout << "Trying factorization of primes below " << upper << " using wheel factorization." << std::endl;

    T candidate = 7;
    uint32_t gapword = WHEEL30GAPS;

    for (size_t i = 0; i < MaxIterations && candidate <= upper; ++i)
    {
        if (primetools::divides(N, candidate)) {
            return std::make_pair(candidate, N / candidate);
        }

        // candidate += gap & 0xf;
        primetools::increment(candidate, gapword & 0xf);

        // Rotate the gaps
        gapword = std::rotr(gapword, 4); // rotate by 4 bits
    }

    return std::nullopt;
}

template <typename T>
const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel210(
    const T& N,
    const size_t MaxIterations
)
{
    if (N < 2) {
        return std::nullopt;
    }
    
    // Get the square root of N as this is our upper bound for trial division
    const T upper = sqrt(N);

    std::cout << "Trying factorization of primes below " << upper << " using wheel factorization." << std::endl;

    T candidate = 7;

    for (size_t i = 0; i < MaxIterations && candidate <= upper; i += (16*3))
    {
        for (auto gapword : WHEEL210GAPS)
        {
            for (size_t index = 0; index < 16; index++)
            {
                if (primetools::divides(N, candidate)) {
                    return std::make_pair(candidate, N / candidate);
                }
                
                // candidate += gap & 0xf;
                primetools::increment(candidate, gapword & 0xf);

                // Rotate the gaps
                gapword = std::rotr(gapword, 4); // rotate by 4 bits
            }
        }
    }

    return std::nullopt;
}

template <typename T>
const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel2310(
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
    while (candidate % 2310 != 1 && candidate > 1) {
        candidate -= 1;
    }

    std::cout << "Trying factorization of " << bits << "-bit primes below " << upper_bound <<
        " using wheel2310 factorization." << std::endl << "Starting candidate: " << candidate << std::endl;

    for (size_t i = 0; i < MaxIterations && candidate <= upper_bound; i += (16*30))
    {
        for (auto gapword : WHEEL2310GAPS)
        {
            for (size_t index = 0; index < 16; index++)
            {
                if (primetools::divides(N, candidate)) {
                    return std::make_pair(candidate, N / candidate);
                }

                // candidate += gap & 0xf;
                primetools::increment(candidate, gapword & 0xf);

                // Rotate the gaps
                gapword = std::rotr(gapword, 4); // rotate by 4 bits
            }
        }
    }

    return std::nullopt;
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
    return TrialDivisionWheel2310<T>(N, MaxIterations);
}

} // namespace primetools

#endif // TRIAL_DIVISION_HPP