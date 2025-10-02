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