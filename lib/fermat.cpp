#include "fermat.hpp"

#include <algorithm>
#include <array>
#include <cinttypes>
#include <iostream>
#include <map>

#include <gmpxx.h>

#include "factors.hpp"
#include "util.hpp"

namespace primetools {

// Given an N to factor, calculate the maximum number of iterations
// to perform in Fermat's factorization method. The given formula is
// The 4th root of N, rounded up to the nearest integer.
const size_t
CalculateFermatIterations(
    const mpz_class& N
)
{
    mpz_class root4;
    mpz_root(root4.get_mpz_t(), N.get_mpz_t(), 4);
    return static_cast<size_t>(root4.get_ui()) + 1;
}

const std::optional<PrimeFactors<mpz_class>>
FermatFactorisation(
    const mpz_class& N,
    const size_t Offset,
    const size_t Max
)
{
    if (N < 2) {
        return std::nullopt;
    }

    mpz_class a, b;

    a = sqrt(N) + 1 + Offset;

    for (size_t i = 0; i < Max; ++i) {
        b = (a * a) - N;

        if (mpz_perfect_square_p(b.get_mpz_t())) {
            return PrimeFactors<mpz_class>::FromPair(
                a - sqrt(b), a + sqrt(b)
            );
        }

        a++;
    }

    return std::nullopt;
}

const std::optional<PrimeFactors<mpz_class>>
FermatFactorisationAlgorithm2(
    const mpz_class& N,
    const size_t Max
)
{
    if (N < 2) {
        return std::nullopt;
    }

    mpz_class u, v, r;
    u = 2 * sqrt(N);
    v = 0;

    for (size_t i = 0; i < Max; ++i) {
        r = (u * u) - (v * v) - 4 * N;
        if (r == 0) {
            mpz_class p = (u - v) / 2;
            mpz_class q = (u + v) / 2;
            return PrimeFactors<mpz_class>::FromPair(p, q);
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
const std::optional<PrimeFactors<mpz_class>>
ModifiedFermatFactorisation4(
    const mpz_class& N,
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

    mpz_class x = sqrt(N) + 1;
    uint64_t NMod20 = mpz_fdiv_ui(N.get_mpz_t(), 20);
    
    // Repeat until x is congruent to 0 mod 10
    // or we found a factor (ChangeN function)
    while (mpz_fdiv_ui(x.get_mpz_t(), 10) != 0) {
        mpz_class b = (x * x) - N;
        if (mpz_perfect_square_p(b.get_mpz_t())) {
            return PrimeFactors<mpz_class>::FromPair(
                x - sqrt(b), x + sqrt(b)
            );
        }
        x++;
    }

    mpz_class y = sqrt(x * x - N);

    uint32_t check = MFFV4XMod10LUT[NMod20];

    for (size_t i = 0; i < Max; i++) {
        const uint64_t x_mod_10 = mpz_fdiv_ui(x.get_mpz_t(), 10);
        if ((1 << x_mod_10) & check) {
            y = x * x - N;
            if (mpz_perfect_square_p(y.get_mpz_t())) {
                mpz_class p = x - sqrt(y);
                mpz_class q = x + sqrt(y);
                return PrimeFactors<mpz_class>::FromPair(p, q);
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
const std::optional<PrimeFactors<mpz_class>>
FMMod20Precomp(
    const mpz_class& N,
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
    while (mpz_fdiv_ui(u.get_mpz_t(), 10) != 0) {
        mpz_class b = (u * u) - N;
        if (mpz_perfect_square_p(b.get_mpz_t())) {
            return PrimeFactors<mpz_class>::FromPair(
                u - sqrt(b), u + sqrt(b)
            );
        }
        u++;
    }

    const uint64_t r = mpz_fdiv_ui(N.get_mpz_t(), 20);
    uint32_t gaps = DIF_LUT.at(r);

    // Implement the cycle routine
    for (size_t i = 0; i < Max; i++)
    {
        primetools::increment(u, gaps & 0xf);
        const mpz_class b = (u * u) - N;
        if (mpz_perfect_square_p(b.get_mpz_t())) {
            return PrimeFactors<mpz_class>::FromPair(
                u - sqrt(b), u + sqrt(b)
            );
        }
        // Rotate the gaps
        gaps = std::rotr(gaps, 4);
    }

    return std::nullopt;
}
    

}