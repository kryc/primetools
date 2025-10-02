#include <array>
#include <bit>
#include <cassert>
#include <iostream>
#include <optional>
#include <future>
#include <vector>
#include <atomic>
#include <thread>

#include "trial_division.hpp"

namespace primetools {

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

    const size_t bits = Bits != 0 ? Bits : (GuessSize ? primetools::GuessSizeOfPrimeFactors(N, true) : primetools::bit_size(N));
    mpz_class lower_bound = GuessSize ? (mpz_class(1) << (bits - 1)) + Range.first : Range.first;
    mpz_class upper_bound = (Range.second == 0) 
        ? (GuessSize ? mpz_class(1) << bits : mpz_class(sqrt(N))) 
        : (GuessSize ? lower_bound + Range.second : Range.second);

    // Step back lower_bound to be a multiple of Modulus
    if (lower_bound % Modulus != 0) {
        lower_bound -= (lower_bound % Modulus);
    }

    // Increment upper_bound to be a multiple of Modulus
    if (upper_bound % Modulus != 0) {
        upper_bound += (Modulus - (upper_bound % Modulus));
    }

    std::atomic<bool> found{false};
    std::vector<std::future<std::optional<std::pair<mpz_class, mpz_class>>>> futures;

    const size_t BlockSize = Modulus * 512;

    for (size_t i = 0; i < num_threads; i++) {
        futures.emplace_back(std::async(std::launch::async, [&, i]() {
            // Each thread runs its own random search
            mpz_class thread_start = (i * BlockSize);
            mpz_class thread_end = (i * BlockSize) + BlockSize;
            while (!found.load(std::memory_order_relaxed)) {
                auto res = TrialDivision<mpz_class, true>(N, Modulus, false, bits, thread_start, thread_end);
                if (res) {
                    found.store(true, std::memory_order_relaxed);
                    return res;
                }
                // Optionally, add a small sleep to avoid busy-waiting
                std::this_thread::yield();
                primetools::increment(thread_start, Threads * BlockSize);
                primetools::increment(thread_end, Threads * BlockSize);
            }
            return std::optional<std::pair<mpz_class, mpz_class>>{};
        }));
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

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionRandomMT(
    const mpz_class& N,
    const size_t Threads,
    const bool GuessSize,
    const size_t Bits,
    const Range<mpz_class>& Range
)
{
    // Use hardware concurrency if Threads == 0
    const size_t num_threads = Threads ? Threads : std::thread::hardware_concurrency();
    const size_t iterations = 1 << 20; // Each thread will process 2^20 iterations before checking for completion

    std::atomic<bool> found{false};
    std::vector<std::future<std::optional<std::pair<mpz_class, mpz_class>>>> futures;

    for (size_t i = 0; i < num_threads; i++) {
        uint64_t seed = (i + 1);
        futures.emplace_back(std::async(std::launch::async, [&, seed]() {
            // Each thread runs its own random search
            auto thread_seed = seed;
            while (!found.load(std::memory_order_relaxed)) {
                auto res = TrialDivisionRandom<5, 510510, 7680, WHEEL510510GAPS, true>(N, GuessSize, Bits, Range, thread_seed, iterations);
                if (res) {
                    found.store(true, std::memory_order_relaxed);
                    return res;
                }
                // Optionally, add a small sleep to avoid busy-waiting
                std::this_thread::yield();
            }
            thread_seed += num_threads; // Change seed for next iteration
            return std::optional<std::pair<mpz_class, mpz_class>>{};
        }));
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

} // namespace primetools