#include <array>
#include <bit>
#include <cassert>
#include <iostream>
#include <numeric>
#include <optional>
#include <future>
#include <vector>
#include <atomic>
#include <thread>
#include <mutex>

#include "prime.hpp"
#include "trial_division.hpp"

namespace primetools {

// Mutex for thread-safe status output
std::mutex status_mutex;

// Worker function for TrialDivisionMT
std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionMTWorker(
    const mpz_class& N,
    const mpz_class& LowerBound,
    const mpz_class& UpperBound,
    const mpz_class& ChunkSize,
    const mpz_class& Chunks,
    const size_t ThreadId,
    const size_t NumThreads,
    const size_t Modulus,
    std::atomic<bool>& Found
) {
    mpz_class thread_start = LowerBound + (ThreadId * ChunkSize);
    mpz_class thread_end = thread_start + ChunkSize;
    while (!Found.load(std::memory_order_relaxed) && thread_start <= UpperBound) {
        mpz_class chunk_index = (thread_start - LowerBound) / ChunkSize;
        {
            std::lock_guard<std::mutex> lock(status_mutex);
            std::cout << '\r' << "Chunk " << chunk_index << " of " << Chunks << " (" <<
                (chunk_index * 100) / Chunks << "%) " << std::flush;
        }
        auto result = TrialDivisionRange<mpz_class>(N, thread_start, thread_end, Modulus);
        if (result) {
            Found.store(true, std::memory_order_relaxed);
            return result;
        }
        std::this_thread::yield();
        primetools::increment(thread_start, NumThreads * ChunkSize);
        primetools::increment(thread_end, NumThreads * ChunkSize);
    }
    return std::nullopt;
}

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionMT(
    const mpz_class& N,
    const size_t Threads,
    const size_t BlockSize,
    const bool GuessSize,
    const size_t Bits,
    const mpz_class& RangeLower,
    const mpz_class& RangeUpper,
    const size_t Modulus
)
{
    // Use hardware concurrency if Threads == 0
    const size_t num_threads = Threads ? Threads : std::thread::hardware_concurrency();
    const size_t block_size = BlockSize ? BlockSize : 512;

    // Get upper and lower bounds
    auto [lower_bound, upper_bound, bits] = GetUpperAndLowerBounds<mpz_class>(N, Modulus, GuessSize, Bits, RangeLower, RangeUpper);

    std::atomic<bool> found{false};
    std::vector<std::future<std::optional<std::pair<mpz_class, mpz_class>>>> futures;

    const mpz_class chunk_size = Modulus * block_size;
    const mpz_class diff = upper_bound - lower_bound;
    const mpz_class chunks = (diff / chunk_size);

    std::cout << "Trying factorization of primes in range [" << primetools::TruncateNumber(lower_bound) << ", " << primetools::TruncateNumber(upper_bound) <<
            "] using modulus " << Modulus << ". " << chunks << " chunks" << std::endl;

    for (size_t i = 0; i < num_threads; i++) {
        futures.emplace_back(std::async(std::launch::async, TrialDivisionMTWorker,
            std::cref(N),
            std::cref(lower_bound),
            std::cref(upper_bound),
            std::cref(chunk_size),
            std::cref(chunks),
            i,
            num_threads,
            Modulus,
            std::ref(found)
        ));
    }

    std::optional<std::pair<mpz_class, mpz_class>> result;
    for (auto& fut : futures) {
        auto res = fut.get();
        if (res) {
            result = res;
            break;
        }
    }

    // Terminate the status line
    std::cout << std::endl;
    return result;
}

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

// Worker function for TrialDivisionRandomMT
std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionRandomMTWorker(
    const mpz_class& N,
    const mpz_class& LowerBound,
    const mpz_class& Chunks,
    const size_t Modulus,
    const size_t ChunkSize,
    const size_t ThreadID,
    std::atomic<bool>& Found
) {
    // Initialize the PRNG
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, ThreadID + 1);

    mpz_class start_block;
    while (!Found.load(std::memory_order_relaxed)) {
        mpz_urandomm(start_block.get_mpz_t(), state, Chunks.get_mpz_t());
        mpz_class thread_start = LowerBound + (start_block * ChunkSize);
        mpz_class thread_end = thread_start + ChunkSize;
        // {
        //     std::lock_guard<std::mutex> lock(status_mutex);
        //     std::cout << '\r' << "Trying block index " << start_block << std::flush;
        // }
        auto result = TrialDivisionRange<mpz_class>(N, thread_start, thread_end, Modulus);
        if (result) {
            Found.store(true, std::memory_order_relaxed);
            return result;
        }
        std::this_thread::yield();
    }
    return std::nullopt;
}

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionRandomMT(
    const mpz_class& N,
    const size_t Threads,
    const size_t BlockSize,
    const bool GuessSize,
    const size_t Bits,
    const mpz_class& RangeLower,
    const mpz_class& RangeUpper,
    const size_t Modulus
)
{
    // Use hardware concurrency if Threads == 0
    const size_t num_threads = Threads ? Threads : std::thread::hardware_concurrency();
    const size_t block_size = BlockSize ? BlockSize : 512;

    std::atomic<bool> found{false};
    std::vector<std::future<std::optional<std::pair<mpz_class, mpz_class>>>> futures;

    // Get bounds and bits
    auto [lower_bound, upper_bound, bits] = GetUpperAndLowerBounds<mpz_class>(N, Modulus, GuessSize, Bits, RangeLower, RangeUpper);

    // Calculate the number of Modulus-sized blocks
    const size_t chunk_size = Modulus * block_size;
    mpz_class diff = upper_bound - lower_bound;
    mpz_class chunks = diff / chunk_size;

    std::cout << "Trying factorization of primes in range [" << primetools::TruncateNumber(lower_bound) << ", " << primetools::TruncateNumber(upper_bound) <<
            "] using random block search with modulus " << Modulus << ". " << chunks << " chunks" << std::endl;

    for (size_t i = 0; i < num_threads; i++) {
        futures.emplace_back(std::async(std::launch::async, TrialDivisionRandomMTWorker,
            std::cref(N),
            std::cref(lower_bound),
            std::cref(chunks),
            Modulus,
            chunk_size,
            i,
            std::ref(found)
        ));
    }

    std::optional<std::pair<mpz_class, mpz_class>> result;
    for (auto& fut : futures) {
        auto res = fut.get();
        if (res) {
            result = res;
            break;
        }
    }
    
    return result;
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
    for (size_t residue = 3; residue < modulus; residue += 2) {
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

} // namespace primetools