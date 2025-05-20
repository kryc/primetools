#ifndef FACTORISE_HPP
#define FACTORISE_HPP

#include <optional>
#include <utility>

#include <gmpxx.h>

#include "fermat.hpp"
#include "pollards_rho.hpp"

namespace primetools {

const std::optional<std::pair<mpz_class, mpz_class>>
FactoriseSmallPrimes(
    const mpz_class& N
);

const std::optional<std::pair<mpz_class, mpz_class>>
Factorise(
    const mpz_class& N
);

}

#endif // FACTORISE_HPP