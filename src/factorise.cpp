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

const std::optional<std::pair<mpz_class, mpz_class>>
FactorisePerfectSquare(
    const mpz_class& N
)
{
    if (N < 2) {
        return std::nullopt;
    }

    if (mpz_perfect_square_p(N.get_mpz_t())) {
        mpz_class sqrtN = sqrt(N);
        return std::make_pair(sqrtN, sqrtN);
    }

    return std::nullopt;
}

// Factorise against the list of the first 100,000 primes
const std::optional<std::pair<mpz_class, mpz_class>>
FactoriseSmallPrimes(
    const mpz_class& N
)
{
    if (N < 2) {
        return std::nullopt;
    }

    unsigned long int prime = 0;

    for (const auto prime_difference : g_FirstPrimes) {
        prime += prime_difference;
        if (mpz_divisible_ui_p(N.get_mpz_t(), prime)) {
            mpz_class factor = N / prime;
            return std::make_pair(prime, factor);
        }
    }

    return std::nullopt;
}

// Factorise against next _Count_ primes
const std::optional<std::pair<mpz_class, mpz_class>>
FactorisePrimesInRange(
    const mpz_class& N,
    const mpz_class& Start,
    const size_t Count
)
{
    if (N < 2) {
        return std::nullopt;
    }

    mpz_class prime = Start;

    for (size_t i = 0; i < Count; ++i) {
        if (N % prime == 0) {
            return std::make_pair(prime, N / prime);
        }
        mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
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

    // Check for perfect square
    std::cout << "Checking for perfect square..." << std::endl;
    auto result = FactorisePerfectSquare(N);
    if (result) {
        return result;
    }

    // Check for small prime factors
    std::cout << "Checking small primes..." << std::endl;
    result = FactoriseSmallPrimes(N);
    if (result) {
        return result;
    }

    // Use Fermat's factorization method up to 2^16 iterations
    std::cout << "Trying Fermat's factorization..." << std::endl;
    result = FermatFactorisation(N, 0, (size_t)1 << 24);
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