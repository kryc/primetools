#ifndef FERMAT_HPP
#define FERMAT_HPP

#include <optional>
#include <utility>

#include <gmpxx.h>

namespace primetools {

const std::optional<std::pair<mpz_class, mpz_class>>
FermatFactorization(
    const mpz_class& N,
    const size_t Max
);

inline const std::optional<std::pair<mpz_class, mpz_class>>
FermatFactorization(
    const mpz_class& N
)
{
    return FermatFactorization(N, (std::numeric_limits<size_t>::max)());
}

}

#endif // FERMAT_HPP