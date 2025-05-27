#include <iostream>
#include <limits>
#include <optional>
#include <span>
#include <utility>

#include <gmpxx.h>

#include "factorise.hpp"
#include "fermat.hpp"
#include "pollards_rho.hpp"

namespace primetools {

static const unsigned char g_Primes[] = {
#embed "../rsrc/small_primes.bin"
};

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

    // primes_data is a vlq-encoded array of the differences between primes
    for (size_t i = 0; i < sizeof(g_Primes);) {
        // Decode vlq-encoded prime difference
        size_t primediff = 0;
        size_t shift = 0;
        uint8_t byte;
        do {
            byte = g_Primes[i++];
            primediff |= (byte & 0x7F) << shift;
            shift += 7;
        } while (i < sizeof(g_Primes) && (byte & 0x80) != 0);

        prime += primediff;

        if (mpz_divisible_ui_p(N.get_mpz_t(), prime)) {
            mpz_class factor = N / prime;
            return std::make_pair(mpz_class(prime), factor);
        }
    }

    return std::nullopt;
}

// Factorise against next _Count_ primes
const std::optional<std::pair<mpz_class, mpz_class>>
FactorisePrimesInRange(
    const mpz_class& N,
    const mpz_class& Start,
    const mpz_class& End
)
{
    if (N < 2) {
        return std::nullopt;
    }

    mpz_class prime;

    // Check if start is prime
    if (mpz_probab_prime_p(prime.get_mpz_t(), 25) == 0) {
        mpz_nextprime(prime.get_mpz_t(), prime.get_mpz_t());
    }

    while (prime <= End) {
        if (mpz_divisible_p(N.get_mpz_t(), prime.get_mpz_t())) {
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

    // Use Fermat's factorization method up to 2^24 iterations
    size_t iterations = CalculateFermatIterations(N);
    iterations = std::min(iterations, (size_t)1 << 24);
    std::cout << "Trying " << iterations << " iterations of Fermat's factorization..." << std::endl;
    result = FermatFactorisation(N, 0, iterations);
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