#ifndef SHANKS_HPP
#define SHANKS_HPP

#include <functional>
#include <optional>
#include <utility>

namespace primetools {

std::optional<std::pair<mpz_class, mpz_class>>
SQUFOF(
    const mpz_class& N,
    const size_t Max = std::numeric_limits<size_t>::max()
);

}

#endif // SHANKS_HPP