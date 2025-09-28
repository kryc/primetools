#include <array>
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
    std::array<uint64_t, 3> wheel210_nibbles = {
        0x462664246264242aULL, // gaps[0..15]
        0x4264684242486462ULL, // gaps[16..31]
        0x2a24246264246626ULL  // gaps[32..47]
    };

    mpz_class candidate = 7;

    for (size_t i = 0; i < MaxIterations && candidate <= upper; i += (16*3))
    {
        for (auto gapword : wheel210_nibbles)
        {
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

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionWheel2310(
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

    // 2310‑wheel gaps packed 16×4‑bit per 64‑bit word
    std::array<uint64_t, 30> wheel2310_nibbles = {
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

    mpz_class candidate = 1;

    for (size_t i = 0; i < MaxIterations && candidate <= upper; i += (16*30))
    {
        for (auto gapword : wheel2310_nibbles)
        {
            for (size_t index = 0; index < 16; index++)
            {
                // Check if it divides N
                if (mpz_divisible_p(N.get_mpz_t(), candidate.get_mpz_t()) && candidate > 1) {
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