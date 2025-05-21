#include <iostream>
#include <limits>
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

    for (const auto& prime : g_FirstPrimes) {
        if (N % prime == 0) {
            return std::make_pair(prime, N / prime);
        }
    }

    return std::nullopt;
}

// Factorise against next _Count_ primes
const std::optional<std::pair<mpz_class, mpz_class>>
FactoriseSmallPrimes2(
    const mpz_class& N,
    const size_t Count
)
{
    if (N < 2) {
        return std::nullopt;
    }

    mpz_class prime = g_FirstPrimes.back();

    for (size_t i = 0; i < Count; ++i) {
        mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
        if (N % prime == 0) {
            return std::make_pair(prime, N / prime);
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
    std::cout << "Checking small primes..." << std::endl;
    auto result = FactoriseSmallPrimes(N);
    if (result) {
        return result;
    }

    // Use Fermat's factorization method up to 2^16 iterations
    std::cout << "Trying Fermat's factorization..." << std::endl;
    result = FermatFactorisation(N, 0, (size_t)1 << 16);
    if (result) {
        return result;
    }

    // Continue small primes
    std::cout << "Checking more small primes..." << std::endl;
    result = FactoriseSmallPrimes2(N, 900000);
    if (result) {
        return result;
    }

    // Try some more Fermat's factorization
    std::cout << "Trying Fermat's factorization again..." << std::endl;
    result = FermatFactorisation(N, (size_t)1 << 16, (size_t)1 << 28);
    if (result) {
        return result;
    }

    // Use Pollard's rho algorithm
    std::cout << "Trying Pollard's rho..." << std::endl;
    result = BrentPollardsRho(N, std::numeric_limits<size_t>::max());
    if (result) {
        return result;
    }

    return std::nullopt;
}

}