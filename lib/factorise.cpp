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
#include "trial_division.hpp"

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

}