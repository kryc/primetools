#ifndef POLLARDS_RHO_HPP
#define POLLARDS_RHO_HPP

#include <algorithm>
#include <atomic>
#include <functional>
#include <mutex>
#include <optional>
#include <thread>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include "factors.hpp"
#include "prime.hpp"
#include "primegenerator.hpp"
#include "util.hpp"

namespace primetools {

namespace {
    // The default max iterations is 2^32
    static const size_t DefaultMaxIterations = (size_t)(1) << 32;
    // The default starting value is 2
    static const mpz_class DefaultStartingValue = 2;
    // The deafult M is 64
    static const size_t DefaultM = 64;

    static const size_t PMinus1DefaultB = (size_t)(1) << 22;
    static const size_t PMinus1DefaultBases = 16ULL;
}

template <typename T>
const T
PollardsRhoPolynomial1(
    const T& x,
    const T& N
)
{
    return (x * x - 1) % N;
}

template <typename T>
const T
PollardsRhoPolynomial2(
    const T& x,
    const T& N
)
{
    return (x * x + 1) % N;
}

template <typename T>
std::optional<std::pair<T, T>>
PollardsRho(
    const T N,
    const std::function<T(T, T)> Polynomial = PollardsRhoPolynomial1<T>,
    const T StartingValue = DefaultStartingValue,
    const size_t Max = DefaultMaxIterations
)
{
    T x = StartingValue;
    T y = x;
    T d = 1;

    for (size_t i = 0; i < Max; ++i) {
        x = Polynomial(x, N);
        y = Polynomial(Polynomial(y, N), N);
        T z = (x > y) ? (x - y) : (y - x); //Avoid the abs call
        // z = primetools::abs(z);
        d = primetools::Gcd(z, N);

        if (d > 1 && d < N) {
            return std::make_pair(d, N / d);
        }
    }

    return std::nullopt;
}

template <typename T>
std::optional<std::pair<T, T>>
BrentPollardsRho(
    const T n,
    const size_t block_size = DefaultM,
    const T starting_value = DefaultStartingValue,
    // NOTE: This is a cap on polynomial evaluations ("rho steps"), not on Brent cycle-doublings.
    // Doubling-based iteration limits grow the work exponentially and are easy to misuse.
    const size_t max_steps = DefaultMaxIterations
)
{
    T current = starting_value;
    T cycle_point = current;
    T gcd_result = 1;
    T product = 1;
    T saved_point = current;
    T diff;

    size_t cycle_length = 1;
    size_t steps = 0;

    const auto f = [&n](const T& x) -> T {
        return (x * x - 1) % n;
    };

    while (steps < max_steps && gcd_result == 1) {
        cycle_point = current;

        for (size_t j = 0; j < cycle_length && steps < max_steps; ++j) {
            current = f(current);
            ++steps;
        }

        size_t block_counter = 0;
        product = 1;

        while (block_counter < cycle_length && gcd_result == 1 && steps < max_steps) {
            saved_point = current;

            const size_t this_block = std::min(block_size, cycle_length - block_counter);
            for (size_t m = 0; m < this_block && gcd_result == 1 && steps < max_steps; ++m) {
                current = f(current);
                ++steps;
                diff = (current > cycle_point) ? (current - cycle_point) : (cycle_point - current); //Avoid the abs call
                // diff = primetools::abs(diff);
                // mpz_abs(diff.get_mpz_t(), diff.get_mpz_t());
                product = (product * diff) % n;
            }
            gcd_result = primetools::Gcd(product, n);
            // mpz_gcd(gcd_result.get_mpz_t(), product.get_mpz_t(), n.get_mpz_t());
            block_counter += this_block;
        }
        cycle_length *= 2;
    }

    if (gcd_result == n) {
        // Standard Brent fallback: walk forward from the saved point until a non-trivial gcd appears.
        // Use gcd(|x - ys|, n) rather than reusing the accumulated product.
        do {
            if (steps >= max_steps) {
                break;
            }
            saved_point = f(saved_point);
            ++steps;

            diff = (cycle_point > saved_point) ? (cycle_point - saved_point) : (saved_point - cycle_point); //Avoid the abs call
            gcd_result = primetools::Gcd(diff, n);
        } while (gcd_result == 1);
    }

    if (gcd_result > 1 && gcd_result < n) {
        return std::make_pair(gcd_result, n / gcd_result);
    }

    return std::nullopt;
}

template <typename T>
std::optional<std::pair<T, T>>
BrentPollardsRhoMT(
    const T N,
    const size_t Threads,
    const size_t M = DefaultM,
    const size_t MaxSteps = DefaultMaxIterations
)
{
    std::vector<std::thread> thread_pool;
    std::atomic<bool> found(false);
    std::optional<std::pair<T, T>> result;
    std::mutex result_mutex;

    const size_t num_threads = Threads ? Threads : std::max<size_t>(1, std::thread::hardware_concurrency());

    // NOTE: PrimeGenerator is not thread-safe (it mutates internal state).
    // Precompute independent starting values on the main thread.
    std::vector<T> starting_values;
    starting_values.reserve(num_threads);
    {
        PrimeGenerator<T> primegen;
        for (size_t i = 0; i < num_threads; ++i) {
            starting_values.emplace_back(T(primegen.Next()));
        }
    }

    for (size_t thread_id = 0; thread_id < num_threads; ++thread_id) {
        const T thread_starting_value = starting_values[thread_id];
        thread_pool.emplace_back([&, thread_starting_value]() {
            auto thread_result = BrentPollardsRho<T>(N, M, thread_starting_value, MaxSteps);
            if (thread_result && !found.exchange(true)) {
                std::lock_guard<std::mutex> lock(result_mutex);
                result = thread_result;
            }
        });
    }

    for (auto& thread : thread_pool) {
        if (thread.joinable()) {
            thread.join();
        }
    }

    return result;
}

template <typename T>
static inline std::optional<std::pair<T, T>>
PollardsPMinus1Stage1(
    const T& N,
    const size_t Bound,
    const size_t Bases,
    const std::span<const uint64_t> Primes,
    const size_t BaseOffset = 0,
    const size_t BaseStride = 1,
    std::atomic<bool>* Found = nullptr
)
{
    for (size_t base_index = BaseOffset; base_index < Bases; base_index += BaseStride) {
        if (Found && Found->load(std::memory_order_relaxed)) {
            return std::nullopt;
        }

        // The base "a" does not need to be prime; it only needs to be coprime to n.
        T a = T(2) + base_index;

        T d = primetools::Gcd<T>(a, N);
        if (d > 1 && d < N) {
            return std::make_pair(d, N / d);
        }

        for (const size_t p : Primes) {
            // Compute p^k <= Bound (largest power of p not exceeding Bound)
            size_t exp = p;
            while (exp <= (Bound / p)) {
                exp *= p;
            }
            a = primetools::ModExp(a, mpz_class(exp), N);
            d = primetools::Gcd<T>(a - 1, N);
            if (d > 1 && d < N) {
                return std::make_pair(d, N / d);
            }
        }
    }

    return std::nullopt;
}

/*
 * Pollard's NP-1 Algorithm
 * This is a special case factorization algorithm that
 * is effective when one of the factors consists only of
 * the product of small primes.
 */
template <typename T>
std::optional<std::pair<T, T>>
PollardsPMinus1(
    const T& N,
    const size_t B = PMinus1DefaultB,
    const size_t Bases = PMinus1DefaultBases
)
{
    if (N < 2) {
        return std::nullopt;
    }
    const auto primes = GetPrimesTo(B);
    return PollardsPMinus1Stage1<T>(N, B, Bases, primes);

    return std::nullopt;
}

template <typename T>
std::optional<std::pair<T, T>>
PollardsPMinus1MT(
    const T& N,
    const size_t Threads = 0,
    const size_t B = PMinus1DefaultB,
    const size_t Bases = PMinus1DefaultBases
)
{
    std::vector<std::thread> thread_pool;
    std::atomic<bool> found(false);
    std::optional<std::pair<T, T>> result;
    std::mutex result_mutex;
    const size_t num_threads = Threads ? Threads : std::thread::hardware_concurrency();
    const auto primes = GetPrimesTo(B);

    for (size_t thread_id = 0; thread_id < num_threads; ++thread_id) {
        thread_pool.emplace_back([thread_id, num_threads, &result, &found, &result_mutex, &N, B, Bases, &primes]() {
            auto thread_result = PollardsPMinus1Stage1<T>(
                N,
                B,
                Bases,
                primes,
                /*BaseOffset=*/thread_id,
                /*BaseStride=*/num_threads,
                /*Found=*/&found
            );

            if (thread_result) {
                std::lock_guard<std::mutex> lock(result_mutex);
                if (!found.exchange(true)) {
                    result = thread_result;
                }
            }
        });
    }

    for (auto& thread : thread_pool) {
        if (thread.joinable()) {
            thread.join();
        }
    }

    return result;
}

} // namespace primetools

#endif // POLLARDS_RHO_HPP