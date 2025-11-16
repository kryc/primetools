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

#include "prime.hpp"
#include "trial_division.hpp"

namespace primetools {

// Mutex for getting generated modulus wheels
std::mutex wheels_mutex;
std::map<size_t, std::vector<__uint128_t>> generated_wheels;

const std::optional<std::pair<mpz_class, mpz_class>>
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

const std::vector<__uint128_t>
GenerateWheelGapsForModulus(
    const size_t Modulus,
    const size_t BitSize,
    const PackingType Packing
)
{
    const size_t WordBits = sizeof(__uint128_t) * 8;
    const size_t GapsPerWord = Packing == PackingType::FastPack ? (WordBits - 1) / BitSize : WordBits / BitSize;

    // Calculate the modulus using the first primes
    // to ensure it is valid
    size_t modulus = 1;
    for (size_t prime = 1; modulus < Modulus; prime++) {
        modulus *= primetools::GetNthPrime(prime);
    }

    if (modulus != Modulus) {
        std::cerr << "Error: Unsupported modulus for trial division: " << Modulus << std::endl;
        return {};
    }

    // Generate the gaps and pack them into uint128_t words
    std::vector<__uint128_t> gaps;
    __uint128_t next = 0;
    size_t count = 0;
    uint64_t last_residue = 1;
    uint32_t wheel30Gaps = std::rotr(WHEEL30GAPSUINT32, 4);
    for (size_t residue = 7; residue < modulus; residue += (wheel30Gaps & 0xf), wheel30Gaps = std::rotr(wheel30Gaps, 4)) {
        if (std::gcd(residue, modulus) != 1) {
            continue;
        }
        // std::cout << residue - last_residue << std::endl;
        const __uint128_t gap = Packing == PackingType::Unpacked ? (residue - last_residue) : (residue - last_residue) >> 1;
        if (gap > (1ULL << BitSize)) {
            std::cerr << "Error: Gap exceeds bit size " << BitSize << " for modulus " << Modulus << std::endl;
            return {};
        }
        const size_t shift = count++ * BitSize;
        next |= gap << shift;
        if (count == GapsPerWord) {
            next <<= Packing == PackingType::FastPack ? 1 : 0;
            gaps.push_back(next);
            next = 0;
            count = 0;
        }
        last_residue = residue;
    }

    // Add the gap from the last residue back to the modulus
    const __uint128_t gap = Packing == PackingType::Unpacked ? (1 + (modulus - last_residue)) : (1 + (modulus - last_residue)) >> 1;
    if (count > 0) {
        const size_t shift = count * BitSize;
        next |= gap << shift;
    }
    else {
        next = gap;
    }

    next <<= Packing == PackingType::FastPack ? 1 : 0;
    gaps.push_back(next);

    return gaps;
}

const std::span<const __uint128_t>
GetWheelGapsForModulus(
    const size_t Modulus
)
{
    std::lock_guard<std::mutex> lock(wheels_mutex);
    auto it = generated_wheels.find(Modulus);
    if (it != generated_wheels.end()) {
        return it->second;
    }
    std::cout << "Generating wheel gaps for modulus " << Modulus << std::endl;
    auto gaps = GenerateWheelGapsForModulus(Modulus, 5, PackingType::FastPack);
    generated_wheels[Modulus] = std::move(gaps);
    return generated_wheels[Modulus];
}

} // namespace primetools