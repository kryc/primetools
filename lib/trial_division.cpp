#include <array>
#include <atomic>
#include <bit>
#include <cassert>
#include <future>
#include <iostream>
#include <map>
#include <mutex>
#include <numeric>
#include <optional>
#include <span>
#include <string>
#include <thread>
#include <vector>

#include <gmpxx.h>

#include "factors.hpp"
#include "prime.hpp"
#include "trial_division.hpp"

namespace primetools {

const std::optional<PrimeFactors<mpz_class>>
TrialDivisionRandom(
    const mpz_class& N,
    const bool GuessSize,
    const size_t Bits,
    const mpz_class& RangeLower,
    const mpz_class& RangeUpper,
    const uint64_t Seed,
    const size_t Modulus,
    const size_t MaxIterations
)
{
    if (N < 2) {
        return std::nullopt;
    }

    // Calculate the bounds and bits
    auto [lower_bound, upper_bound, bits] = GetUpperAndLowerBounds<mpz_class>(N, Modulus, GuessSize, Bits, RangeLower, RangeUpper);

    const size_t ChunkSize = Modulus * 512;
    // Calculate the difference between the bounds
    mpz_class diff = upper_bound - lower_bound;

    // Calculate the number of blocks
    mpz_class chunks = diff / ChunkSize;

    std::cout << "Trying random factorization of primes of size " << bits << " bits. " << chunks << " chunks." << std::endl;
    std::cout << "WARNING: This is highly inefficient and not recommended!" << std::endl;

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, Seed);

    mpz_class block;
    mpz_urandomm(block.get_mpz_t(), state, chunks.get_mpz_t());

    for (size_t i = 0; i < MaxIterations; i++) {
        mpz_class block_start = lower_bound + (block * ChunkSize);
        mpz_class block_end = block_start + ChunkSize;

        auto result = TrialDivisionRange<mpz_class>(N, block_start, block_end, Modulus);
        if (result) {
            return result;
        }
        
        // Move to a new random block
        mpz_urandomm(block.get_mpz_t(), state, chunks.get_mpz_t());
    }

    return std::nullopt;
}

} // namespace primetools