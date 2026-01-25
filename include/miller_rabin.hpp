#ifndef MILLERRABIN_HPP
#define MILLERRABIN_HPP

#include <cstdint>

#include "random.hpp"
#include "util.hpp"

#include <gmpxx.h>

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