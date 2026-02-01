#ifndef FERMAT_HPP
#define FERMAT_HPP

#include <array>
#include <atomic>
#include <bit>
#include <cinttypes>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <optional>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include "factors.hpp"

namespace primetools {

constexpr size_t FermatDefaultOffset = 0;
constexpr size_t FermatDefaultCount = std::numeric_limits<size_t>::max();
constexpr size_t FermatDefaultStride = 1;

typedef enum FermatAlgorithm {
    AlgFermat,
    AlgFermat2,
    AlgModifiedFermatV4,
    AlgFMMod20Precomp
} FermatAlgorithm;

const FermatAlgorithm
GetFermatAlgorithmFromString(
    const std::string_view AlgorithmStr
);

const std::string_view
FermatAlgorithmToString(
    const FermatAlgorithm Algorithm
);

// Given an N to factor, calculate the maximum number of iterations
// to perform in Fermat's factorization method. The given formula is
// The 4th root of N, rounded up to the nearest integer.
template <typename T>
const size_t
CalculateFermatIterations(
    const T& N
)
{
    if constexpr (std::is_same_v<T, mpz_class>) {
        mpz_class root4;
        mpz_root(root4.get_mpz_t(), N.get_mpz_t(), 4);
        return static_cast<size_t>(root4.get_ui()) + 1;
    } else {
        T root4 = primetools::Root(N, 4);
        return static_cast<size_t>(root4) + 1;
    }
}

template <typename T>
const std::optional<std::pair<T, T>>
FermatFactorisationAlgorithm1(
    const T& N,
    const size_t Offset,
    const size_t Max = std::numeric_limits<size_t>::max()
)
{
    if (N < 2) {
        return std::nullopt;
    }

    T a, b;

    a = primetools::Sqrt(N) + 1 + Offset;

    for (size_t i = 0; i < Max; ++i) {
        b = (a * a) - N;

        if (primetools::IsPerfectSquare(b)) {
            return std::make_pair(
                a - primetools::Sqrt(b),
                a + primetools::Sqrt(b)
            );
        }

        a++;
    }

    return std::nullopt;
}

template <typename T>
const std::optional<std::pair<T, T>>
FermatFactorisationAlgorithm2(
    const T& N,
    const size_t Max
)
{
    if (N < 2) {
        return std::nullopt;
    }

    T u, v, r;
    u = 2 * primetools::Sqrt(N);
    v = 0;

    for (size_t i = 0; i < Max; ++i) {
        r = (u * u) - (v * v) - 4 * N;
        if (r == 0) {
            T p = (u - v) / 2;
            T q = (u + v) / 2;
            return std::make_pair(p, q);
        }
        else if (r > 0) {
            v += 2;
        } else {
            u += 2;
        }
    }

    return std::nullopt;
}

/*
 * Modified Fermat Factorisation Version 4 (MFFV4)
 * This version is based on Somsuk's MFFV4 algorithm as described in
 * "A New Modified Integer Factorization Algorithm using Integer Modulo 20's Technique"
 * https://ieeexplore.ieee.org/document/6978214
 */
template <typename T>
const std::optional<std::pair<T, T>>
ModifiedFermatFactorisation4(
    const T& N,
    const size_t Max
)
{
    if (N < 2) {
        return std::nullopt;
    }

    static const std::array<uint32_t,20> MFFV4XMod10LUT = {
        0,       0b1000100010,  // 1 = { 1, 5, 9 }
        0,       0b0100000100,  // 3 = { 2, 8 }
        0, 0, 0, 0b0001010000,  // 7 = { 4, 6 }
        0,       0b0010101000,  // 9 = { 3, 5, 7 }
        0,       0b0001010001,  // 11 = { 0, 4, 6 }
        0,       0b0010001000,  // 13 = { 3, 7 }
        0, 0, 0, 0b1000000010,  // 17 = { 1, 9 }
        0,       0b0100000101   // 19 = { 0, 2, 8 }
    };

    T x = primetools::Sqrt(N) + 1;
    uint64_t NMod20 = primetools::modulo(N, 20);
    
    // Repeat until x is congruent to 0 mod 10
    // or we found a factor (ChangeN function)
    while (primetools::modulo(x, 10) != 0) {
        T b = (x * x) - N;
        if (primetools::IsPerfectSquare(b)) {
            return std::make_pair(
                x - primetools::Sqrt(b), x + primetools::Sqrt(b)
            );
        }
        x++;
    }

    mpz_class y = primetools::Sqrt(x * x - N);

    uint32_t check = MFFV4XMod10LUT[NMod20];

    for (size_t i = 0; i < Max; i++) {
        const uint64_t x_mod_10 = primetools::modulo(x, 10);
        if ((1 << x_mod_10) & check) {
            y = x * x - N;
            if (primetools::IsPerfectSquare(y)) {
                T p = x - primetools::Sqrt(y);
                T q = x + primetools::Sqrt(y);
                return std::make_pair(p, q);
            }
        }
        x++;
    }

    return std::nullopt;
}

/*
 * FMMod20Precomp
 * This algorithm is based on Hatem M. Bahig's paper
 * "Speeding Up Fermat’s Factoring Method using Precomputation"
 * It works based on the observation that
 * The value of 𝑢 𝑚𝑜𝑑 10 (or LSD(𝑢)) is variable, but for a fixed
 * value of 𝑛 𝑚𝑜𝑑 20, the possible values of u, PS(𝑢2 − 𝑛) = True,
 * can be determined initially.
 */
template <typename T>
const std::optional<std::pair<T, T>>
FMMod20Precomp(
    const T& N,
    const size_t Count = FermatDefaultCount,
    const size_t Offset = FermatDefaultOffset,
    const size_t Stride = FermatDefaultStride
)
{
    if (N < 2) {
        return std::nullopt;
    }

    // This algorithm assumes the factors 2 and 5 have been removed
    // Handle that gracefully here
    if (primetools::modulo(N, 2) == 0) {
        return std::make_pair(T(2), N / T(2));
    }
    if (primetools::modulo(N, 5) == 0) {
        return std::make_pair(T(5), N / T(5));
    }

    // Increment gaps. Zero gaps are padding and should be skipped.
    static const std::map<uint64_t, uint64_t> DIF_LUT = {
        {1,  0x1441144114411441}, //  1 = { 1, 4, 4, 1 }
        {3,  0x2062206220622062}, //  3 = { 2, 6, 0, 2 }
        {7,  0x4024402440244024}, //  7 = { 4, 2, 0, 4 }
        {9,  0x3223322332233223}, //  9 = { 3, 2, 2, 3 }
        {11, 0x4240424042404240}, // 11 = { 0, 4, 2, 4 }
        {13, 0x3043304330433043}, // 13 = { 3, 4, 0, 3 }
        {17, 0x1081108110811081}, // 17 = { 1, 8, 0, 1 }
        {19, 0x2620262026202620}  // 19 = { 0, 2, 6, 2 }
    };

    constexpr size_t kBitsPerGap = 4; // Fixed
    constexpr uint64_t kGapMask = (1 << kBitsPerGap) - 1;
    constexpr size_t kValuesPerCycle = sizeof(uint64_t) * 8 / kBitsPerGap; // Fixed (16)
    constexpr size_t kCyclesPerBlock = 64; // Can be changed to tune. Must be a multitple of kValuesPerCycle (16).
    constexpr size_t kValuesPerBlock = kCyclesPerBlock * kValuesPerCycle;
    constexpr size_t kIncrementPerCycle = 40; // Each cycle increases u by 40
    constexpr size_t kIncrementPerBlock = kIncrementPerCycle * kCyclesPerBlock;

    T u0 = primetools::Sqrt(N) + 1;

    // Repeat until u is congruent to 0 mod 10
    // (or we found a factor)
    while (primetools::modulo(u0, 10) != 0) {
        const T b = (u0 * u0) - N;
        if (primetools::IsPerfectSquare(b)) {
            return std::make_pair(
                u0 - primetools::Sqrt(b), u0 + primetools::Sqrt(b)
            );
        }
        u0++;
    }

    const uint64_t r = primetools::modulo(N, 20);
    uint64_t gaps = DIF_LUT.at(r);

    // Apply offset/stride in whole blocks (each block advances u by kIncrementPerBlock).
    // This partitions cleanly for multithreading: thread t uses Offset=t, Stride=num_threads.
    T u = u0 + T(Offset) * T(kIncrementPerBlock);

    // Check the thread's starting point as well (important when Offset > 0).
    {
        const T b = (u * u) - N;
        if (primetools::IsPerfectSquare(b)) {
            return std::make_pair(
                u - primetools::Sqrt(b), u + primetools::Sqrt(b)
            );
        }
    }

    for (size_t block = 0; block < Count; block++)
    {
        // Each block is kCyclesPerBlock cycles, where each cycle consumes kValuesPerCycle gaps.
        for (size_t i = 0; i < kValuesPerBlock; ++i) {
            const uint64_t gap = gaps & kGapMask;
            gaps = std::rotr(gaps, kBitsPerGap);
            if (gap == 0) {
                continue;
            }

            primetools::increment(u, gap);

            const T b = (u * u) - N;
            if (primetools::IsPerfectSquare(b)) {
                return std::make_pair(
                    u - primetools::Sqrt(b), u + primetools::Sqrt(b)
                );
            }
        }

        // Skip (Step-1) whole blocks to avoid overlap with other threads.
        if (Stride > 1) {
            primetools::increment(u, T(kIncrementPerBlock) * T(Stride - 1));
        }
    }

    return std::nullopt;
}

template <typename T>
const std::optional<std::pair<T, T>>
FMMod20PrecompMT(
    const T& N,
    const size_t Threads = 0,
    const size_t Count = FermatDefaultCount,
    const size_t Offset = FermatDefaultOffset
)
{
    std::vector<std::thread> thread_pool;
    std::atomic<bool> found(false);
    std::optional<std::pair<T, T>> result;
    std::mutex result_mutex;

    const size_t num_threads = Threads ? Threads : std::thread::hardware_concurrency();

    for (size_t thread_id = 0; thread_id < num_threads; ++thread_id) {
        thread_pool.emplace_back([thread_id, num_threads, &result, &found, &result_mutex, &N, Count, Offset]() {
            auto thread_result = FMMod20Precomp<T>(
                N,
                Count,
                Offset + thread_id,
                num_threads /* stride */
            );
            if (thread_result && !found.load()) {
                std::lock_guard<std::mutex> lock(result_mutex);
                if (!found.load()) {
                    found.store(true);
                    result = thread_result;
                }
            }
        });
    }

    for (auto& thread : thread_pool) {
        thread.join();
    }

    return result;
}

// Public interface to factorise a number N
template <typename T>
const std::optional<PrimeFactors<T>>
FermatFactorisation(
    const T& N,
    const size_t Threads = 0,
    const FermatAlgorithm Algorithm = AlgFMMod20Precomp,
    const size_t Offset = FermatDefaultOffset,
    const size_t Count = FermatDefaultCount
)
{
    const size_t num_threads = Threads ? Threads : std::thread::hardware_concurrency();
    
    switch (Algorithm) {
        case AlgFermat: {
            auto result = FermatFactorisationAlgorithm1(N, Offset, Count);
            if (result) {
                return PrimeFactors<T>::FromPair(result.value());
            }
            return std::nullopt;
        }
        case AlgFermat2: {
            auto result = FermatFactorisationAlgorithm2(N, Count);
            if (result) {
                return PrimeFactors<T>::FromPair(result->first, result->second);
            }
            return std::nullopt;
        }
        case AlgModifiedFermatV4: {
            auto result = ModifiedFermatFactorisation4(N, Count);
            if (result) {
                return PrimeFactors<T>::FromPair(result->first, result->second);
            }
            return std::nullopt;
        }
        case AlgFMMod20Precomp: {
            std::optional<std::pair<T, T>> result;
            result = num_threads > 1 ? FMMod20PrecompMT(N, num_threads, Count, Offset) : FMMod20Precomp(N, Count, Offset, 1);
            if (result) {
                return PrimeFactors<T>::FromPair(result->first, result->second);
            }
            return std::nullopt;
        }
        default:
            return std::nullopt;
    }
}

} // namespace primetools

#endif // FERMAT_HPP