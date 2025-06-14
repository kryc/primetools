#include <iostream>
#include <limits>
#include <optional>
#include <random>
#include <span>
#include <utility>

#include <gmpxx.h>

#include "analyse.hpp"
#include "factorise.hpp"
#include "fermat.hpp"
#include "pollard.hpp"
#include "random.hpp"

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
    std::cout << "Trying " << iterations << " iterations of FMMod20Precomp (Fermat) factorization..." << std::endl;
    result = FMMod20Precomp(N, iterations);
    if (result) {
        return result;
    }

    // Try using Pollards P-1
    std::cout << "Trying Pollard's P-1 (B2**20)..." << std::endl;
    result = PollardsPMinus1(N, (size_t)1 << 20);
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

const std::optional<std::pair<mpz_class, mpz_class>>
RandomPrimeFactorization(
    const mpz_class& N,
    const size_t Base,
    const size_t MaxIterations
)
{
    if (N < 2) {
        return std::nullopt;
    }

    // Calculate the size of N's constituent primes
    const size_t bits = primetools::GuessSizeOfPrimeFactors(N, true);

    std::cout << "Trying random factorization of primes of size " << bits << " bits." << std::endl;
    std::cout << "WARNING: This is highly inefficient and not recommended!" << std::endl;
    
    // Calculate the lower and upper bounds for the prime factor
    const mpz_class lower_bound = mpz_class(1) << (bits - 1);
    const mpz_class upper_bound = (mpz_class(1) << bits);

    // Calculate the difference between the bounds
    const mpz_class diff = upper_bound - lower_bound;

    // We now have a base offset that will bounce around the
    // range of lower and upper bounds.
    mpz_class base = Base;

    // Set up our PRNG
    primetools::MiniPRNG64 prng;

    for (size_t i = 0; i < MaxIterations; ++i) {
        // Generate a random candidate prime
        const mpz_class candidate = lower_bound + base;

        // Check if it divides N
        if (mpz_divisible_p(N.get_mpz_t(), candidate.get_mpz_t())) {
            return std::make_pair(candidate, N / candidate);
        }

        // Pseudorandomly adjust the base (2a + 1) % diff.
        // The choice of (2a + 1) is arbitrary but ensures
        // that base remains an odd.
        base = (base + prng.NextEven()) % diff;
    }

    return std::nullopt;
}

}