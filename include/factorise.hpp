#ifndef FACTORISE_HPP
#define FACTORISE_HPP

#include <optional>
#include <utility>

#include <gmpxx.h>

#include "factors.hpp"
#include "fermat.hpp"
#include "pollard.hpp"
#include "shanks.hpp"
#include "trial_division.hpp"

namespace primetools {

// Factorise a perfect square
const std::optional<PrimeFactors<mpz_class>>
FactorisePerfectSquare(
    const mpz_class& N
);

// Public interface to factorise a number N
template <typename T>
const std::optional<PrimeFactors<T>>
Fermat(
    const T& N,
    const FermatAlgorithm Algorithm,
    const size_t Offset,
    const size_t Max = std::numeric_limits<size_t>::max()
)
{
    switch (Algorithm) {
        case AlgFermat: {
            auto result = FermatFactorisation(N, Offset, Max);
            if (result) {
                return PrimeFactors<T>::FromPair(result.value());
            }
            return std::nullopt;
        }
        case AlgFermat2: {
            auto result = FermatFactorisationAlgorithm2(N, Max);
            if (result) {
                return PrimeFactors<T>::FromPair(result->first, result->second);
            }
            return std::nullopt;
        }
        case AlgModifiedFermatV4: {
            auto result = ModifiedFermatFactorisation4(N, Max);
            if (result) {
                return PrimeFactors<T>::FromPair(result->first, result->second);
            }
            return std::nullopt;
        }
        case AlgFMMod20Precomp: {
            auto result = FMMod20Precomp(N, Max);
            if (result) {
                return PrimeFactors<T>::FromPair(result->first, result->second);
            }
            return std::nullopt;
        }
        default:
            return std::nullopt;
    }
}

static inline void Log(
    std::string_view Message,
    const bool Verbose
)
{
    if (Verbose)
        std::cout << Message << std::endl;
}

// Try to factorise with increasing complexity
template <typename T>
const std::optional<PrimeFactors<T>>
Factorise(
    const T& N,
    const size_t Threads = 0,
    const bool Verbose = false
)
{
    if (N < 2) {
        return std::nullopt;
    }

    PrimeFactors<T> factors;
    size_t threads = Threads > 0 ? Threads : std::thread::hardware_concurrency();

    // Check for perfect square
    Log("Checking for perfect square...", Verbose);
    auto result = FactorisePerfectSquare(N);
    if (result) {
        return result;
    }

    // Check for small prime factors
    Log("Checking small primes...", Verbose);
    size_t modulus = 9699690;
    result = TrialDivision(N, 0, 0, false, 0, T(0), T(modulus * threads), modulus);
    if (result) {
        factors = result.value();
    }

    T remainder = N / result->Product();
    if (remainder == 1) {
        return factors;
    }

    Log("Current factors: " + factors.GetString(), Verbose);
    Log("Remainder after small primes: " + remainder.get_str(), Verbose);

    // Use Fermat's factorization method up to 2^24 iterations
    size_t iterations = CalculateFermatIterations(N);
    iterations = std::min(iterations, (size_t)(1) << 24);
    Log("Trying " + std::to_string(iterations) + " iterations of FMMod20Precomp (Fermat) factorization...", Verbose);
    auto resultpair = FMMod20Precomp(remainder, iterations);
    if (resultpair) {
        // If either factor is prime, add it to the factors
        if (primetools::isprime(resultpair->first)) {
            factors.AddFactor(resultpair->first);
        }
        if (primetools::isprime(resultpair->second)) {
            factors.AddFactor(resultpair->second);
        }
    }

    remainder = N / factors.Product();
    if (remainder == 1) {
        return factors;
    }

    Log("Current factors: " + factors.GetString(), Verbose);
    Log("Remainder after FMMod20Precomp: " + remainder.get_str(), Verbose);

    // Try using Pollards P-1
    Log("Trying Pollard's P-1 (B2**32)...", Verbose);
    resultpair = PollardsPMinus1(remainder, (size_t)1 << 22, 100);
    if (resultpair) {
        // If either factor is prime, add it to the factors
        if (primetools::isprime(resultpair->first)) {
            factors.AddFactor(resultpair->first);
        }
        if (primetools::isprime(resultpair->second)) {
            factors.AddFactor(resultpair->second);
        }
    }

    remainder = N / factors.Product();
    if (remainder == 1) {
        return factors;
    }

    Log("Current factors: " + factors.GetString(), Verbose);
    Log("Remainder after P-1: " + remainder.get_str(), Verbose);

    // Use Pollard's rho algorithm
    Log("Trying Brent-Pollard's rho...", Verbose);
    resultpair = BrentPollardsRho(remainder);
    if (resultpair) {
        // If either factor is prime, add it to the factors
        if (primetools::isprime(resultpair->first)) {
            factors.AddFactor(resultpair->first);
        }
        if (primetools::isprime(resultpair->second)) {
            factors.AddFactor(resultpair->second);
        }
    }

    // Fall back to trial division if still not fully factored
    remainder = N / factors.Product();
    if (remainder == 1) {
        return factors;
    }

    Log("Current factors: " + factors.GetString(), Verbose);
    Log("Remainder after Brent-Pollard's rho: " + remainder.get_str(), Verbose);

    Log("Falling back to trial division...", Verbose);
    result = TrialDivision(remainder, threads, 0, false, 0, T(modulus * threads), T(0), modulus);
    if (result) {
        factors.Update(result.value());
        return factors;
    }

    return std::nullopt;
}

}

#endif // FACTORISE_HPP