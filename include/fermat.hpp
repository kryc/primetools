#ifndef FERMAT_HPP
#define FERMAT_HPP

#include <array>
#include <cinttypes>
#include <iostream>
#include <limits>
#include <optional>
#include <string_view>
#include <utility>

#include <gmpxx.h>

#include "factors.hpp"

namespace primetools {

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

const size_t
CalculateFermatIterations(
    const mpz_class& N
);

template <typename T>
const std::optional<std::pair<T, T>>
FermatFactorisation(
    const T& N,
    const size_t Offset,
    const size_t Max = std::numeric_limits<size_t>::max()
)
{
    if (N < 2) {
        return std::nullopt;
    }

    T a, b;

    a = sqrt(N) + 1 + Offset;

    for (size_t i = 0; i < Max; ++i) {
        b = (a * a) - N;

        if (primetools::IsPerfectSquare(b)) {
            return std::make_pair(
                a - sqrt(b),
                a + sqrt(b)
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
    u = 2 * sqrt(N);
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

    T x = sqrt(N) + 1;
    uint64_t NMod20 = primetools::modulo(N, 20);
    
    // Repeat until x is congruent to 0 mod 10
    // or we found a factor (ChangeN function)
    while (primetools::modulo(x, 10) != 0) {
        T b = (x * x) - N;
        if (primetools::IsPerfectSquare(b)) {
            return std::make_pair(
                x - sqrt(b), x + sqrt(b)
            );
        }
        x++;
    }

    mpz_class y = sqrt(x * x - N);

    uint32_t check = MFFV4XMod10LUT[NMod20];

    for (size_t i = 0; i < Max; i++) {
        const uint64_t x_mod_10 = primetools::modulo(x, 10);
        if ((1 << x_mod_10) & check) {
            y = x * x - N;
            if (primetools::IsPerfectSquare(y)) {
                T p = x - sqrt(y);
                T q = x + sqrt(y);
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
    const size_t Max
)
{
    if (N < 2) {
        return std::nullopt;
    }

    static const std::map<uint64_t, uint64_t> DIF_LUT = {
        {1,  0x1441144114411441}, //  1 = { 1, 4, 4, 1 }
        {3,  0x2062206220622062}, //  3 = { 2, 6, -1, 2 }
        {7,  0x4024402440244024}, //  7 = { 4, 2, -1, 4 }
        {9,  0x3223322332233223}, //  9 = { 3, 2, 2, 3 }
        {11, 0x4240424042404240}, // 11 = { 0, 4, 2, 4 }
        {13, 0x3043304330433043}, // 13 = { 3, 4, -1, 3 }
        {17, 0x1081108110811081}, // 17 = { 1, 8, -1, 1 }
        {19, 0x2620262026202620}  // 19 = { 0, 2, 6, 2 }
    };

    mpz_class u = sqrt(N) + 1;

    // Repeat until u is congruent to 0 mod 10
    // (or we found a factor)
    while (primetools::modulo(u, 10) != 0) {
        T b = (u * u) - N;
        if (primetools::IsPerfectSquare(b)) {
            return std::make_pair(
                u - sqrt(b), u + sqrt(b)
            );
        }
        u++;
    }

    const uint64_t r = primetools::modulo(N, 20);
    uint32_t gaps = DIF_LUT.at(r);

    // Implement the cycle routine
    for (size_t i = 0; i < Max; i++)
    {
        primetools::increment(u, gaps & 0xf);
        const T b = (u * u) - N;
        if (primetools::IsPerfectSquare(b)) {
            return std::make_pair(
                u - sqrt(b), u + sqrt(b)
            );
        }
        // Rotate the gaps
        gaps = std::rotr(gaps, 4);
    }

    return std::nullopt;
}

} // namespace primetools

#endif // FERMAT_HPP