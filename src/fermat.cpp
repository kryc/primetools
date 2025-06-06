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

const mpz_class
ChangeN(
    const mpz_class& X,
    const unsigned char NMod20
)
{
    mpz_class x = X;
    switch(NMod20)
    {
    case 1:
        while (x % 10 != 1 && x % 10 != 5 && x % 10 != 9)
            x++;
        break;
    case 3:
        while (x % 10 != 2 && x % 10 != 8)
            x++;
        break;
    case 7:
        while (x % 10 != 4 && x % 10 != 6)
            x++;
        break;
    case 9:
        while (x % 10 != 3 && x % 10 != 5 && x % 10 != 7)
            x++;
        break;
    case 11:
        while (x % 10 != 0 && x % 10 != 4 && x % 10 != 6)
            x++;
        break;
    case 13:
        while (x % 10 != 3 && x % 10 != 7)
            x++;
        break;
    case 17:
        while (x % 10 != 1 && x % 10 != 9)
            x++;
        break;
    case 19:
        while (x % 10 != 0 && x % 10 != 2 && x % 10 != 8)
            x++;
        break;
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

    static const int8_t MFFV4_LUT[20][3] = {
        {},         { 1, 5, 9 },  // 1
        {},         { 2, 8, -1 }, // 3
        {}, {}, {}, { 4, 6, -1 }, // 7
        {}, {}, {}, { 0, 4, -1 }, // 11
        {},         { 3, 7, -1 }, // 13
        {}, {}, {}, { 1, 9, -1 }, // 17
        {},         { 0, 2, 8 }   // 19
    };

    mpz_class x = sqrt(N) + 1;
    unsigned long int NMod20 = mpz_fdiv_ui(N.get_mpz_t(), 20);
    x = ChangeN(x, NMod20);
    mpz_class y = sqrt(x * x - N);

    int8_t check[3];
    check[0] = MFFV4_LUT[NMod20][0];
    check[1] = MFFV4_LUT[NMod20][1];
    check[2] = MFFV4_LUT[NMod20][2];

    for (size_t i = 0; i < Max; i++) {
        const uint64_t x_mod_10 = mpz_fdiv_ui(x.get_mpz_t(), 10);
        if (x_mod_10 == check[0] ||
            x_mod_10 == check[1] ||
            x_mod_10 == check[2]) {
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

}