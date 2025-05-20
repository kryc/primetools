#include <optional>
#include <utility>

#include <gmpxx.h>

#include "factorise.hpp"
#include "fermat.hpp"
#include "pollards_rho.hpp"
#include "primes.hpp"

namespace primetools {

// Factorise against the list of the first 100,000 primes
const std::optional<std::pair<mpz_class, mpz_class>>
FactoriseSmallPrimes(
    const mpz_class& N
)
{
    if (N < 2) {
        return std::nullopt;
    }

    mpz_class n = N;
    mpz_class factor = 0;

    for (const auto& prime : g_FirstPrimes) {
        if (n % prime == 0) {
            return std::make_pair(prime, n / prime);
        }
    }

    return std::nullopt;
}

const std::optional<std::pair<mpz_class, mpz_class>>
Factorise(
    const mpz_class& N
)
{
    if (N < 2) {
        return std::nullopt;
    }

    // Check for small prime factors
    auto result = FactoriseSmallPrimes(N);
    if (result) {
        return result;
    }

    // Use Fermat's factorization method up to 2^20 iterations
    result = FermatFactorization(N, 1 << 20);
    if (result) {
        return result;
    }

    // Use Pollard's rho algorithm
    result = BrentPollardsRho(N);
    if (result) {
        return result;
    }

    return std::nullopt;
}

}