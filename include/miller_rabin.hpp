#ifndef MILLERRABIN_HPP
#define MILLERRABIN_HPP

#include <cstdint>

#include "random.hpp"
#include "util.hpp"

#include <gmpxx.h>

static constexpr size_t kMillerRabinSmallThreshold = 1122004669633;

template <typename T>
const bool
MillerRabinSmall(
    const T& N
)
{
    static constexpr std::array<uint64_t, 2> kWitnessesSub1373653 = {2, 3};
    static constexpr std::array<uint64_t, 2> kWitnessesSub9080191 = {31, 73};
    static constexpr std::array<uint64_t, 3> kWitnessesSub25326001 = {2, 3, 5};
    static constexpr std::array<uint64_t, 7> kWitnessesSub3215031751 = {2, 3, 5, 7, 11, 13, 17};
    static constexpr std::array<uint64_t, 4> kWitnessesSub1122004669633 = {2, 13, 23, 1662803};
    std::span<const uint64_t> witnesses;

    if (N < 1373653) {
        witnesses = kWitnessesSub1373653;
    } else if (N < 9080191) {
        witnesses = kWitnessesSub9080191;
    } else if (N < 25326001) {
        witnesses = kWitnessesSub25326001;
    } else if (N < 3215031751) {
        witnesses = kWitnessesSub3215031751;
    } else if (N < 1122004669633) {
        witnesses = kWitnessesSub1122004669633;
    } else {
        throw std::invalid_argument("N is too large for MillerRabinSmall");
    }

    // Find d such that N = 2^r * d + 1 with d odd
    T d = N - 1;
    size_t r = 0;
    while (primetools::IsEven(d)) {
        d /= 2;
        r++;
    }

    // Witness loop
    for (const auto& a_w : witnesses) {
        T a = static_cast<T>(a_w);
        T x = primetools::ModExp(a, d, N);
        if (x == 1 || x == N - 1) {
            continue;
        }

        bool composite = true;
        for (size_t j = 0; j < r - 1; j++) {
            x = primetools::ModExp(x, 2, N);
            if (x == N - 1) {
                composite = false;
                break;
            }
        }
        if (composite) {
            return false;
        }
    }

    return true;
}

template <typename T>
const bool
MillerRabin(
    const T& N,
    const size_t K
)
{
    if (N < 2) {
        return false;
    }
    if (N != 2 && primetools::IsEven(N)) {
        return false;
    }
    if (N <= 3) {
        return true;
    }

    if (N < kMillerRabinSmallThreshold) {
        return MillerRabinSmall(N);
    }

    // Find d such that N = 2^r * d + 1 with d odd
    T d = N - 1;
    size_t r = 0;
    while (primetools::IsEven(d)) {
        d /= 2;
        r++;
    }

    // Instantiate a PRNG
    primetools::MiniPRNG64 prng;
    const T randModulus = N - 4;

    // Witness loop
    for (size_t i = 0; i < K; i++) {
        T a = primetools::RandomInRange(prng, randModulus) + 2; // Random base in [2, N-2]
        T x = primetools::ModExp(a, d, N);
        if (x == 1 || x == N - 1) {
            continue;
        }

        bool composite = true;
        for (size_t j = 0; j < r - 1; j++) {
            x = primetools::ModExp(x, 2, N);
            if (x == N - 1) {
                composite = false;
                break;
            }
        }
        if (composite) {
            return false;
        }
    }

    return true;
}

#endif // MILLERRABIN_HPP