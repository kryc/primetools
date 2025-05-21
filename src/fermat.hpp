#ifndef FERMAT_HPP
#define FERMAT_HPP

#include <optional>
#include <utility>

#include <gmpxx.h>

namespace primetools {

const std::optional<std::pair<mpz_class, mpz_class>>
FermatFactorisation(
    const mpz_class& N,
    const size_t Offset,
    const size_t Max
);

inline const std::optional<std::pair<mpz_class, mpz_class>>
FermatFactorisation(
    const mpz_class& N
)
{
    return FermatFactorisation(
        N,
        0,
        (std::numeric_limits<size_t>::max)()
    );
}

// Sieve improvement for fermat's factorization
const std::optional<std::pair<mpz_class, mpz_class>>
FermatSieve(
    const mpz_class& N,
    const size_t Offset,
    const size_t Max
);

inline const std::optional<std::pair<mpz_class, mpz_class>>
FermatSieve(
    const mpz_class& N
)
{
    return FermatSieve(
        N,
        0,
        (std::numeric_limits<size_t>::max)()
    );
}

}

#endif // FERMAT_HPP