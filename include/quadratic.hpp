#ifndef QUDRATIC_HPP
#define QUDRATIC_HPP

#include <cstddef>
#include <cstdint>
#include <optional>
#include <utility>

#include <gmpxx.h>

#include "logging.hpp"
#include "util.hpp"

namespace primetools {

// Heuristic: choose a smoothness bound (factor base prime limit) based on the size of N.
// Intended for semiprimes.
const uint32_t
QuadraticSieveRecommendedFactorBaseBound(
    const mpz_class& N
);

// Returns a non-trivial factor pair (p, q) such that p*q == N.
// Returns std::nullopt if no factor was found.
const std::optional<std::pair<mpz_class, mpz_class>>
QuadraticSieveFactor(
    const mpz_class& N,
    LogCallback LogFn = LogStdOut
);
static inline
const std::optional<std::pair<mpz_class, mpz_class>>
QuadraticSieveFactor(
    const __uint128_t& N,
    LogCallback LogFn = LogStdOut
) {
    return QuadraticSieveFactor(primetools::ConvertType<mpz_class>(N), LogFn);
}

} // namespace primetools

#endif // QUDRATIC_HPP
