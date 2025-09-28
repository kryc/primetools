#include <bit>
#include <cassert>
#include <iostream>
#include <optional>


#include "trial_division.hpp"

namespace primetools {

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionBitflip(
    const mpz_class& N,
    const size_t MaxIterations
)
{
    if (N < 2) {
        return std::nullopt;
    }

    // Get the square root of N as this is our upper bound for trial division
    const mpz_class upper = sqrt(N);

    // Get the number of bits in the upper bound
    const size_t bits = mpz_sizeinbase(upper.get_mpz_t(), 2);
    const size_t bit_range = bits - 1;

    std::cout << "Trying random factorization of primes below " << upper.get_str() << std::endl;
    std::cout << "WARNING: This is highly inefficient and not recommended!" << std::endl;
    
    // Make the first candidate odd
    mpz_class candidate = upper | 1;

    // Set up our PRNG
    primetools::MiniPRNG32 prng;

    for (size_t i = 0; i < MaxIterations; ++i) {
        // Flip a random bit in the candidate
        const size_t bit = (prng.Next() % bit_range) + 1;
        candidate ^= (mpz_class(1) << bit);

        // Check if it divides N
        if (mpz_divisible_p(N.get_mpz_t(), candidate.get_mpz_t())) {
            return std::make_pair(candidate, N / candidate);
        }
    }

    return std::nullopt;
}

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel30(
    const mpz_class& N,
    const size_t MaxIterations
)
{
    if (N < 2) {
        return std::nullopt;
    }
    
    // Get the square root of N as this is our upper bound for trial division
    const mpz_class upper = sqrt(N);

    std::cout << "Trying factorization of primes below " << upper.get_str() << " using wheel factorization." << std::endl;

    const uint32_t wheel30gaps = 0x26424246; // {6,4,2,4,2,4,6,2}

    mpz_class candidate = 7;
    uint32_t gap = wheel30gaps;

    for (size_t i = 0; i < MaxIterations && candidate <= upper; ++i) {
        // Check if it divides N
        if (mpz_divisible_p(N.get_mpz_t(), candidate.get_mpz_t())) {
            return std::make_pair(candidate, N / candidate);
        }

        // candidate += gap & 0xf;
        mpz_add_ui(candidate.get_mpz_t(), candidate.get_mpz_t(), gap & 0xf);

        // Rotate the gaps
        gap = std::rotr(gap, 4);
    }

    return std::nullopt;
}

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel210(
    const mpz_class& N,
    const size_t MaxIterations
)
{
    if (N < 2) {
        return std::nullopt;
    }
    
    // Get the square root of N as this is our upper bound for trial division
    const mpz_class upper = sqrt(N);

    std::cout << "Trying factorization of primes below " << upper.get_str() << " using wheel factorization." << std::endl;

    // 210‑wheel gaps packed 16×4‑bit per 64‑bit word
    const uint64_t wheel210_nibbles[3] = {
        0x4626624626424AULL, // gaps[0..15]
        0x42862864268642ULL, // gaps[16..31]
        0x2A242642662642ULL  // gaps[32..47]
    };

    mpz_class candidate = 7;

    for (size_t i = 0; i < MaxIterations && candidate <= upper; i += (16*3))
    {
        for (size_t block = 0; block < 3; block++)
        {
            uint64_t gapword = wheel210_nibbles[block];

            for (size_t index = 0; index < 16; index++)
            {
                // Check if it divides N
                if (mpz_divisible_p(N.get_mpz_t(), candidate.get_mpz_t())) {
                    return std::make_pair(candidate, N / candidate);
                }

                // candidate += gap & 0xf;
                mpz_add_ui(candidate.get_mpz_t(), candidate.get_mpz_t(), gapword & 0xf);

                // Rotate the gaps
                gapword = std::rotr(gapword, 4); // rotate by 4 bits
            }
        }
    }

    return std::nullopt;
}

} // namespace primetools