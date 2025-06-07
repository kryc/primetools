#include "fermat.hpp"

#include <algorithm>
#include <array>
#include <cinttypes>
#include <iostream>


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

const std::optional<std::pair<mpz_class, mpz_class>>
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
            return std::make_pair(a - sqrt(b), a + sqrt(b));
        }

        a++;
    }

    return std::nullopt;
}

const std::optional<std::pair<mpz_class, mpz_class>>
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

const mpz_class
ChangeN(
    const mpz_class& X,
    const unsigned char NMod20
)
{
    const uint32_t check = MFFV4XMod10LUT[NMod20];

    mpz_class x = X;
    uint64_t x_mod_10 = mpz_fdiv_ui(x.get_mpz_t(), 10);

    while (!((1 << x_mod_10) & check))
    {
        x++;
        x_mod_10 = mpz_fdiv_ui(x.get_mpz_t(), 10);
    }

    return x;
}

/*
 * Modified Fermat Factorisation Version 4 (MFFV4)
 * This version is based on Somsuk's MFFV4 algorithm as described in
 * "A New Modified Integer Factorization Algorithm using Integer Modulo 20's Technique"
 * https://ieeexplore.ieee.org/document/6978214
 */
const std::optional<std::pair<mpz_class, mpz_class>>
ModifiedFermatFactorisation4(
    const mpz_class& N,
    const size_t Max
)
{
    if (N < 2) {
        return std::nullopt;
    }

    mpz_class x = sqrt(N) + 1;
    uint64_t NMod20 = mpz_fdiv_ui(N.get_mpz_t(), 20);
    x = ChangeN(x, NMod20);
    mpz_class y = sqrt(x * x - N);

    uint32_t check = MFFV4XMod10LUT[NMod20];

    for (size_t i = 0; i < Max; i++) {
        const uint64_t x_mod_10 = mpz_fdiv_ui(x.get_mpz_t(), 10);
        if ((1 << x_mod_10) & check) {
            y = x * x - N;
            if (mpz_perfect_square_p(y.get_mpz_t())) {
                mpz_class p = x - sqrt(y);
                mpz_class q = x + sqrt(y);
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
const std::optional<std::pair<mpz_class, mpz_class>>
FMMod20Precomp(
    const mpz_class& N,
    const size_t Max
)
{
    if (N < 2) {
        return std::nullopt;
    }

    static const uint32_t DIF_LUT[20] = {
        0,       0x01040401, //  1 = { 1, 4, 4, 1 }
        0,       0x02ff0602, //  3 = { 2, 6, -1, 2 }
        0, 0, 0, 0x04ff0204, //  7 = { 4, 2, -1, 4 }
        0,       0x03020203, //  9 = { 3, 2, 2, 3 }
        0,       0x04020400, // 11 = { 0, 4, 2, 4 }
        0,       0x03ff0403, // 13 = { 3, 4, -1, 3 }
        0, 0, 0, 0x01ff0801, // 17 = { 1, 8, -1, 1 }
        0,       0x02060200  // 19 = { 0, 2, 6, 2 }
    };

    mpz_class u = sqrt(N) + 1;

    // Repeat until u is congruent to 0 mod 10
    // (or we found a factor)
    while (mpz_fdiv_ui(u.get_mpz_t(), 10) != 0) {
        mpz_class b = (u * u) - N;
        if (mpz_perfect_square_p(b.get_mpz_t())) {
            return std::make_pair(u - sqrt(b), u + sqrt(b));
        }
        u++;
    }

    const uint64_t r = mpz_fdiv_ui(N.get_mpz_t(), 20);
    const uint32_t dif = DIF_LUT[r];

    for (size_t i = 0; i < Max; i++)
    {
        uint32_t d = dif;
        // Implement the cycle routine
        for (size_t j = 0; j < 4; j++, d >>= 8) {
            if ((d & 0xff) == 0xff) {
                continue;
            }
            u += (d & 0xff);
            const mpz_class b = (u * u) - N;
            if (mpz_perfect_square_p(b.get_mpz_t())) {
                return std::make_pair(u - sqrt(b), u + sqrt(b));
            }
        }
    }

    return std::nullopt;
}
    

}