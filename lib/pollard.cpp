#include <functional>
#include <optional>
#include <utility>

#include "gmpxx.h"

#include "pollard.hpp"
#include "primegenerator.hpp"
#include "util.hpp"

namespace primetools {

const mpz_class
PollardsRhoPolynomial1(
    const mpz_class& x,
    const mpz_class& N
)
{
    return (x * x - 1) % N;
}

const mpz_class
PollardsRhoPolynomial2(
    const mpz_class& x,
    const mpz_class& N
)
{
    return (x * x + 1) % N;
}

std::optional<std::pair<mpz_class, mpz_class>>
PollardsRho(
    const mpz_class N,
    const std::function<mpz_class(mpz_class, mpz_class)> Polynomial,
    const mpz_class StartingValue,
    const size_t Max
)
{
    mpz_class x = StartingValue;
    mpz_class y = x;
    mpz_class d = 1;

    for (size_t i = 0; i < Max; ++i) {
        x = Polynomial(x, N);
        y = Polynomial(Polynomial(y, N), N);
        mpz_class z = x - y;
        z = primetools::abs(z);
        d = primetools::gcd(z, N);

        if (d > 1 && d < N) {
            return std::make_pair(d, N / d);
        }
    }

    return std::nullopt;
}

std::optional<std::pair<mpz_class, mpz_class>>
BrentPollardsRho(
    const mpz_class N,
    const size_t M,
    const mpz_class StartingValue,
    const size_t Max
)
{
    mpz_class x = StartingValue;
    mpz_class y = x;
    mpz_class d = 1;
    mpz_class q = 1;
    mpz_class k, xs, z;
    mpz_class l = 1;

    for (size_t i = 0; i < Max && d == 1; ++i) {
        y = x;
        for (size_t i = 0; i < l; ++i) {
            x = (x * x - 1) % N;
        }
        k = 0;
        while (k < l && d == 1) {
            xs = x;
            for (size_t i = 0; i < M && d == 1; ++i) {
                x = (x * x - 1) % N;
                z = y - x;
                // z = primetools::abs(z);
                mpz_abs(z.get_mpz_t(), z.get_mpz_t());
                q = (q * z) % N;
            }
            // d = primetools::gcd(q, N);
            mpz_gcd(d.get_mpz_t(), q.get_mpz_t(), N.get_mpz_t());
            k += M;
        }
        l *= 2;
    }

    if (d == N) {
        do {
            xs = (xs * xs - 1) % N;
            z = y - xs;
            // z = primetools::abs(z);
            mpz_abs(z.get_mpz_t(), z.get_mpz_t());
            q = (q * z) % N;
            // d = primetools::gcd(q, N);
            mpz_gcd(d.get_mpz_t(), q.get_mpz_t(), N.get_mpz_t());
        } while (d == 1);
    }

    if (d > 1 && d < N) {
        return std::make_pair(d, N / d);
    }

    return std::nullopt;
}

/*
 * Pollard's NP-1 Algorithm
 * This is a special case factorization algorithm that
 * is effective when one of the factors consists only of
 * the product of small primes.
 */
std::optional<std::pair<mpz_class, mpz_class>>
PollardsPMinus1(
    const mpz_class& N,
    const size_t B,
    const size_t Bases
)
{
    if (N < 2) {
        return std::nullopt;
    }

    PrimeGenerator<mpz_class> primegen;
    PrimeGenerator<mpz_class> baseprimegen;

    for (size_t base = 1; base < Bases; ++base) {
        mpz_class a = baseprimegen.Next();
        mpz_class d;
        mpz_class p = primegen.Next();

        while (p <= B) {
            // Compute the highest power of p that is <= B
            mpz_class exp = p;
            while (exp * p <= B) {
                exp *= p;
            }
            a = primetools::modexp(a, exp, N);
            d = primetools::gcd(a - 1, N);
            if (d > 1 && d < N) {
                return std::make_pair(d, N / d);
            }
            p = primegen.Next();
        }
    }

    return std::nullopt;
}

} // namespace primetools