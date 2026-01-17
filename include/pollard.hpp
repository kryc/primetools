#ifndef POLLARDS_RHO_HPP
#define POLLARDS_RHO_HPP

#include <functional>
#include <optional>
#include <utility>

#include <gmpxx.h>

#include "factors.hpp"
#include "primegenerator.hpp"
#include "util.hpp"

namespace primetools {

namespace {
    // The default max iterations is 2^32
    static const size_t DefaultMaxIterations = (size_t)1 << 32;
    // The default starting value is 2
    static const mpz_class DefaultStartingValue = 2;
    // The deafult M is 1000
    static const size_t DefaultM = 1000;
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
    const T N,
    const size_t M = DefaultM,
    const T StartingValue = DefaultStartingValue,
    const size_t Max = std::numeric_limits<size_t>::max()
)
{
    T x = StartingValue;
    T y = x;
    T d = 1;
    T q = 1;
    T k, xs, z;
    T l = 1;

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
                z = (x > y) ? (x - y) : (y - x); //Avoid the abs call
                // z = primetools::abs(z);
                // mpz_abs(z.get_mpz_t(), z.get_mpz_t());
                q = (q * z) % N;
            }
            d = primetools::gcd(q, N);
            // mpz_gcd(d.get_mpz_t(), q.get_mpz_t(), N.get_mpz_t());
            k += M;
        }
        l *= 2;
    }

    if (d == N) {
        do {
            xs = (xs * xs - 1) % N;
            z = (y > xs) ? (y - xs) : (xs - y); //Avoid the abs call
            // z = primetools::abs(z);
            // mpz_abs(z.get_mpz_t(), z.get_mpz_t());
            q = (q * z) % N;
            d = primetools::gcd(q, N);
            // mpz_gcd(d.get_mpz_t(), q.get_mpz_t(), N.get_mpz_t());
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
template <typename T>
std::optional<std::pair<T, T>>
PollardsPMinus1(
    const T& N,
    const size_t B,
    const size_t Bases
)
{
    if (N < 2) {
        return std::nullopt;
    }

    PrimeGenerator<T> primegen;
    PrimeGenerator<T> baseprimegen;

    for (size_t base = 1; base < Bases; ++base) {
        T a = baseprimegen.Next();
        T d;
        T p = primegen.Next();

        while (p <= B) {
            // Compute the highest power of p that is <= B
            T exp = p;
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

}

#endif // POLLARDS_RHO_HPP