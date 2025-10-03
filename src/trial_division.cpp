#include <array>
#include <bit>
#include <cassert>
#include <iostream>
#include <optional>
#include <future>
#include <vector>
#include <atomic>
#include <thread>
#include <mutex>

#include "trial_division.hpp"

namespace primetools {

// Mutex for thread-safe status output
std::mutex status_mutex;

// Worker function for TrialDivisionMT
std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionMTWorker(
    const mpz_class& N,
    const mpz_class& lower_bound,
    const mpz_class& upper_bound,
    const mpz_class& BlockSize,
    const mpz_class& blocks,
    const size_t thread_id,
    const size_t num_threads,
    std::atomic<bool>& found
) {
    mpz_class thread_start = lower_bound + (thread_id * BlockSize);
    mpz_class thread_end = thread_start + BlockSize;
    while (!found.load(std::memory_order_relaxed) && thread_start <= upper_bound) {
        mpz_class block_index = (thread_start - lower_bound) / BlockSize;
        {
            std::lock_guard<std::mutex> lock(status_mutex);
            std::cout << '\r' << "Block " << block_index << " of " << blocks << " (" <<
                (block_index * 100) / blocks << "%) " << std::flush;
        }
        auto res = TrialDivisionRange<mpz_class, 510510, 5, 7680, WHEEL510510GAPS>(N, thread_start, thread_end);
        if (res) {
            found.store(true, std::memory_order_relaxed);
            return res;
        }
        std::this_thread::yield();
        primetools::increment(thread_start, num_threads * BlockSize);
        primetools::increment(thread_end, num_threads * BlockSize);
    }
    return std::nullopt;
}

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionMT(
    const mpz_class& N,
    const size_t Threads,
    const bool GuessSize,
    const size_t Bits,
    const Range<mpz_class>& Range,
    const size_t Modulus
)
{
    // Use hardware concurrency if Threads == 0
    const size_t num_threads = Threads ? Threads : std::thread::hardware_concurrency();

    // Get upper and lower bounds
    auto [lower_bound, upper_bound, bits] = GetUpperAndLowerBounds<mpz_class>(N, Modulus, GuessSize, Bits, Range);

    std::atomic<bool> found{false};
    std::vector<std::future<std::optional<std::pair<mpz_class, mpz_class>>>> futures;

    const mpz_class BlockSize = Modulus * 256;
    const mpz_class diff = upper_bound - lower_bound;
    const mpz_class blocks = (diff / BlockSize);

    std::cout << "Trying factorization of primes in range [" << lower_bound << ", " << upper_bound <<
            "] using modulus " << Modulus << ". " << blocks << " blocks" << std::endl;

    for (size_t i = 0; i < num_threads; i++) {
        futures.emplace_back(std::async(std::launch::async, TrialDivisionMTWorker,
            std::cref(N),
            std::cref(lower_bound),
            std::cref(upper_bound),
            std::cref(BlockSize),
            std::cref(blocks),
            i,
            num_threads,
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

// Worker function for TrialDivisionRandomMT
std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionRandomMTWorker(
    const mpz_class& N,
    const mpz_class& lower_bound,
    const mpz_class& upper_bound,
    const mpz_class& blocks,
    const size_t Modulus,
    const size_t BlockSize,
    const size_t thread_id,
    std::atomic<bool>& found
) {
    // Initialize the PRNG
    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, thread_id + 1);

    mpz_class start_block;
    while (!found.load(std::memory_order_relaxed)) {
        mpz_urandomm(start_block.get_mpz_t(), state, blocks.get_mpz_t());
        mpz_class thread_start = lower_bound + (start_block * BlockSize);
        mpz_class thread_end = thread_start + BlockSize;
        // {
        //     std::lock_guard<std::mutex> lock(status_mutex);
        //     std::cout << '\r' << "Trying block index " << start_block << std::flush;
        // }
        auto res = TrialDivisionRange<mpz_class>(N, thread_start, thread_end, Modulus);
        if (res) {
            found.store(true, std::memory_order_relaxed);
            return res;
        }
        std::this_thread::yield();
    }
    return std::nullopt;
}

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionRandomMT(
    const mpz_class& N,
    const size_t Threads,
    const bool GuessSize,
    const size_t Bits,
    const Range<mpz_class>& Range,
    const size_t Modulus
)
{
    // Use hardware concurrency if Threads == 0
    const size_t BlockSize = Modulus * 2048;
    const size_t num_threads = Threads ? Threads : std::thread::hardware_concurrency();

    std::atomic<bool> found{false};
    std::vector<std::future<std::optional<std::pair<mpz_class, mpz_class>>>> futures;

    // Get bounds and bits
    auto [lower_bound, upper_bound, bits] = GetUpperAndLowerBounds<mpz_class>(N, Modulus, GuessSize, Bits, Range);

    // Calculate the number of Modulus-sized blocks
    mpz_class diff = upper_bound - lower_bound;
    mpz_class blocks = diff / BlockSize;

    std::cout << "Trying factorization of primes in range [" << lower_bound << ", " << upper_bound <<
            "] using random block search with modulus " << Modulus << ". " << blocks << " blocks" << std::endl;

    for (size_t i = 0; i < num_threads; i++) {
        futures.emplace_back(std::async(std::launch::async, TrialDivisionRandomMTWorker,
            std::cref(N),
            std::cref(lower_bound),
            std::cref(upper_bound),
            std::cref(blocks),
            Modulus,
            BlockSize,
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
    // Terminate the status line
    std::cout << std::endl;
    return result;
}

} // namespace primetools