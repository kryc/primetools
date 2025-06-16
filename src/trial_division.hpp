#ifndef TRIAL_DIVISION_HPP
#define TRIAL_DIVISION_HPP

#include <iostream>
#include <optional>

#include <gmpxx.h>

#include "analyse.hpp"
#include "random.hpp"

namespace primetools {

template <bool Random>
const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivision(
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

        if constexpr (Random) {
            base = (base + prng.NextEven()) % diff;
        } else {
            base = (base + 2) % diff;
        }
    }

    return std::nullopt;
}

// Define template specializations for Random and Non-Random cases
inline std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionRandom(const mpz_class& N, const size_t Base, const size_t MaxIterations) {
    return TrialDivision<true>(N, Base, MaxIterations);
}

inline std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionLinear(const mpz_class& N, const size_t Base, const size_t MaxIterations) {
    return TrialDivision<false>(N, Base, MaxIterations);
}

} // namespace primetools

#endif // TRIAL_DIVISION_HPP