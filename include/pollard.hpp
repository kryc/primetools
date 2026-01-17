#ifndef POLLARDS_RHO_HPP
#define POLLARDS_RHO_HPP

#include <atomic>
#include <functional>
#include <mutex>
#include <optional>
#include <thread>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include "factors.hpp"
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
    static const size_t PMinus1DefaultBases = (size_t)(1) << 24;
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
        d = primetools::gcd(z, N);

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
    const size_t max_iterations = DefaultMaxIterations
)
{
    T current = starting_value;
    T cycle_point = current;
    T gcd_result = 1;
    T product = 1;
    T block_counter, saved_point, diff;
    T cycle_length = 1;

    for (size_t iter = 0; iter < max_iterations && gcd_result == 1; ++iter) {
        cycle_point = current;
        for (size_t j = 0; j < cycle_length; ++j) {
            current = (current * current - 1) % n;
        }
        block_counter = 0;
        while (block_counter < cycle_length && gcd_result == 1) {
            saved_point = current;
            for (size_t m = 0; m < block_size && gcd_result == 1; ++m) {
                current = (current * current - 1) % n;
                diff = (current > cycle_point) ? (current - cycle_point) : (cycle_point - current); //Avoid the abs call
                // diff = primetools::abs(diff);
                // mpz_abs(diff.get_mpz_t(), diff.get_mpz_t());
                product = (product * diff) % n;
            }
            gcd_result = primetools::gcd(product, n);
            // mpz_gcd(gcd_result.get_mpz_t(), product.get_mpz_t(), n.get_mpz_t());
            block_counter += block_size;
        }
        cycle_length *= 2;
    }

    if (gcd_result == n) {
        do {
            saved_point = (saved_point * saved_point - 1) % n;
            diff = (cycle_point > saved_point) ? (cycle_point - saved_point) : (saved_point - cycle_point); //Avoid the abs call
            // diff = primetools::abs(diff);
            // mpz_abs(diff.get_mpz_t(), diff.get_mpz_t());
            product = (product * diff) % n;
            gcd_result = primetools::gcd(product, n);
            // mpz_gcd(gcd_result.get_mpz_t(), product.get_mpz_t(), n.get_mpz_t());
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
    const size_t MaxIterations = DefaultMaxIterations
)
{
    std::vector<std::thread> thread_pool;
    std::atomic<bool> found(false);
    std::optional<std::pair<T, T>> result;
    std::mutex result_mutex;
    PrimeGenerator<T> primegen;

    for (size_t thread_id = 0; thread_id < Threads; ++thread_id) {
        thread_pool.emplace_back([&]() {
            T thread_starting_value = primegen.Next();
            auto thread_result = BrentPollardsRho<T>(N, M, thread_starting_value, MaxIterations);
            if (thread_result && !found.load()) {
                found.store(true);
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


// Precompute primes up to the stage-1 bound once.
static inline std::vector<size_t>
PMinus1PrimesUpTo(
    const size_t bound
)
{
    std::vector<size_t> primes;
    primes.reserve(1u << 20);
    PrimeGenerator<size_t> primegen;
    for (size_t p = primegen.Next(); p <= bound; p = primegen.Next()) {
        primes.push_back(p);
    }
    return primes;
}

template <typename T>
static inline std::optional<std::pair<T, T>>
PollardsPMinus1Stage1(
    const T& n,
    const size_t bound,
    const size_t bases,
    const std::vector<size_t>& primes,
    const size_t base_offset = 0,
    const size_t base_stride = 1,
    std::atomic<bool>* stop = nullptr
)
{
    for (size_t base_index = base_offset; base_index < bases; base_index += base_stride) {
        if (stop && stop->load(std::memory_order_relaxed)) {
            return std::nullopt;
        }

        // The base "a" does not need to be prime; it only needs to be coprime to n.
        T a = T(2) + base_index;

        T d = primetools::gcd(a, n);
        if (d > 1 && d < n) {
            return std::make_pair(d, n / d);
        }

        for (const size_t p : primes) {
            // Compute p^k <= bound (largest power of p not exceeding bound)
            size_t exp = p;
            while (exp <= (bound / p)) {
                exp *= p;
            }
            a = primetools::modexp(a, mpz_class(exp), n);
            d = primetools::gcd(a - 1, n);
            if (d > 1 && d < n) {
                return std::make_pair(d, n / d);
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
    const auto primes = PMinus1PrimesUpTo(B);
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
    const auto primes = PMinus1PrimesUpTo(B);

    for (size_t thread_id = 0; thread_id < num_threads; ++thread_id) {
        thread_pool.emplace_back([thread_id, num_threads, &result, &found, &result_mutex, &N, B, Bases, &primes]() {
            auto thread_result = PollardsPMinus1Stage1<T>(
                N,
                B,
                Bases,
                primes,
                /*base_offset=*/thread_id,
                /*base_stride=*/num_threads,
                /*stop=*/&found
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