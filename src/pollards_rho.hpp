#ifndef POLLARDS_RHO_HPP
#define POLLARDS_RHO_HPP

#include <functional>
#include <optional>
#include <utility>

#include "util.hpp"

namespace primetools {

std::optional<std::pair<mpz_class, mpz_class>>
PollardsRho(
    const mpz_class N,
    const std::function<mpz_class(mpz_class, mpz_class)> Polynomial,
    const mpz_class StartingValue,
    const size_t Max
);

std::optional<std::pair<mpz_class, mpz_class>>
PollardsRho(
    const mpz_class N
);

std::optional<std::pair<mpz_class, mpz_class>>
BrentPollardsRho(
    const mpz_class N,
    const size_t M,
    const mpz_class StartingValue,
    const size_t Max
);

std::optional<std::pair<mpz_class, mpz_class>>
BrentPollardsRho(
    const mpz_class N
);

}

#endif // POLLARDS_RHO_HPP