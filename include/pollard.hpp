#ifndef POLLARDS_RHO_HPP
#define POLLARDS_RHO_HPP

#include <functional>
#include <optional>
#include <utility>

#include <gmpxx.h>

#include "factors.hpp"
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

const mpz_class
PollardsRhoPolynomial1(
    const mpz_class& x,
    const mpz_class& N
);

const mpz_class
PollardsRhoPolynomial2(
    const mpz_class& x,
    const mpz_class& N
);

std::optional<PrimeFactors<mpz_class>>
PollardsRho(
    const mpz_class N,
    const std::function<mpz_class(mpz_class, mpz_class)> Polynomial = PollardsRhoPolynomial2,
    const mpz_class StartingValue = DefaultStartingValue,
    const size_t Max = std::numeric_limits<size_t>::max()
);

std::optional<PrimeFactors<mpz_class>>
BrentPollardsRho(
    const mpz_class N,
    const size_t M = DefaultM,
    const mpz_class StartingValue = DefaultStartingValue,
    const size_t Max = std::numeric_limits<size_t>::max()
);

std::optional<PrimeFactors<mpz_class>>
PollardsPMinus1(
    const mpz_class& N,
    const size_t B = (size_t)1 << 20,
    const size_t Bases = 1
);

}

#endif // POLLARDS_RHO_HPP