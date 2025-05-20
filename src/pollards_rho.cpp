#include <functional>
#include <optional>
#include <utility>

#include "gmpxx.h"

#include "pollards_rho.hpp"
#include "util.hpp"

namespace {

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

}

namespace primetools {

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
PollardsRho(
    const mpz_class N
)
{
    // The original algoritm used a polynomial of the form x^2 - 1
    // but more common now is x^2 + 1
    static const auto DefaultPolynomial = std::function<mpz_class(mpz_class, mpz_class)>(PollardsRhoPolynomial2);
    // The default max iterations is 2^32
    static const size_t DefaultMaxIterations = (size_t)1 << 32;
    // The default starting value is 2
    static const mpz_class DefaultStartingValue = 2;
    
    return PollardsRho(
        N,
        DefaultPolynomial,
        DefaultStartingValue,
        DefaultMaxIterations
    );
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
            x = PollardsRhoPolynomial2(x, N);
        }
        k = 0;
        while (k < l && d == 1) {
            xs = x;
            for (size_t i = 0; i < M && d == 1; ++i) {
                x = PollardsRhoPolynomial2(x, N);
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
            xs = PollardsRhoPolynomial2(xs, N);
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

std::optional<std::pair<mpz_class, mpz_class>>
BrentPollardsRho(
    const mpz_class N
)
{
    // The default max iterations is 2^32
    static const size_t DefaultMaxIterations = (size_t)1 << 32;
    // The default starting value is 2
    static const mpz_class DefaultStartingValue = 2;
    // The deafult M is 1000
    static const size_t DefaultM = 1000;
    
    return BrentPollardsRho(
        N,
        DefaultM,
        DefaultStartingValue,
        DefaultMaxIterations
    );
}

} // namespace primetools