#ifndef FACTORISE_HPP
#define FACTORISE_HPP

#include <optional>
#include <utility>

#include <gmpxx.h>

#include "fermat.hpp"
#include "pollard.hpp"
#include "shanks.hpp"
#include "trial_division.hpp"

namespace primetools {

// Factorise a perfect square
const std::optional<std::pair<mpz_class, mpz_class>>
FactorisePerfectSquare(
    const mpz_class& N
);

// Try to factorise with increasing complexity
template <typename T>
const std::optional<std::pair<mpz_class, mpz_class>>
Factorise(
    const T& N
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
    result = TrialDivisionWheel510510(N, 1000000, false);
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
    std::cout << "Trying Pollard's P-1 (B2**22)..." << std::endl;
    result = PollardsPMinus1(N, (size_t)1 << 22);
    if (result) {
        return result;
    }

    // Use Pollard's rho algorithm
    std::cout << "Trying Brent-Pollard's rho..." << std::endl;
    result = BrentPollardsRho(N, std::numeric_limits<size_t>::max());
    if (result) {
        return result;
    }

    return std::nullopt;
}

}

#endif // FACTORISE_HPP