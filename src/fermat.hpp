#ifndef FERMAT_HPP
#define FERMAT_HPP

#include <limits>
#include <optional>
#include <utility>

#include <gmpxx.h>

namespace primetools {

const size_t
CalculateFermatIterations(
    const mpz_class& N
);

const std::optional<std::pair<mpz_class, mpz_class>>
FermatFactorisation(
    const mpz_class& N,
    const size_t Offset = 0,
    const size_t Max = std::numeric_limits<size_t>::max()
);

const std::optional<std::pair<mpz_class, mpz_class>>
FermatFactorisationAlgorithm2(
    const mpz_class& N,
    const size_t Max = std::numeric_limits<size_t>::max()
);

const std::optional<std::pair<mpz_class, mpz_class>>
ModifiedFermatFactorisation4(
    const mpz_class& N,
    const size_t Max = std::numeric_limits<size_t>::max()
);

}

#endif // FERMAT_HPP