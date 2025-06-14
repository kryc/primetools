#ifndef FACTORISE_HPP
#define FACTORISE_HPP

#include <optional>
#include <utility>

#include <gmpxx.h>

#include "fermat.hpp"
#include "pollard.hpp"
#include "shanks.hpp"

namespace primetools {

// Factorise a perfect square
const std::optional<std::pair<mpz_class, mpz_class>>
FactorisePerfectSquare(
    const mpz_class& N
);

// Factorise a number with small prime factors
const std::optional<std::pair<mpz_class, mpz_class>>
FactoriseSmallPrimes(
    const mpz_class& N
);

// Try to factorise with increasing complexity
const std::optional<std::pair<mpz_class, mpz_class>>
Factorise(
    const mpz_class& N
);

// Completely random candidate factor selection
// This is not a recommended method for factorisation
const std::optional<std::pair<mpz_class, mpz_class>>
RandomPrimeFactorization(
    const mpz_class& N,
    const size_t Base = 1,
    const size_t MaxIterations = std::numeric_limits<size_t>::max()
);

}

#endif // FACTORISE_HPP