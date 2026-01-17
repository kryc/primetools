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

// Try to factorise with increasing complexity
template <typename T>
const std::optional<PrimeFactors<T>>
Factorise(
    const T& N,
    const size_t Threads = 0
)
{
    if (N < 2) {
        return std::nullopt;
    }

    PrimeFactors<T> factors;
    size_t threads = Threads > 0 ? Threads : std::thread::hardware_concurrency();

    // Check for perfect square
    std::cout << "Checking for perfect square..." << std::endl;
    auto result = FactorisePerfectSquare(N);
    if (result) {
        return result;
    }

    // Check for small prime factors
    std::cout << "Checking small primes..." << std::endl;
    size_t modulus = 9699690;
    result = TrialDivision(N, 0, 0, false, 0, T(0), T(modulus * threads), modulus);
    if (result) {
        factors = result.value();
    }

    T remainder = N / result->Product();
    if (remainder == 1) {
        return factors;
    }

    std::cout << "Current factors: " << factors.GetString() << std::endl;
    std::cout << "Remainder after small primes: " << remainder << std::endl;

    // Use Fermat's factorization method up to 2^24 iterations
    size_t iterations = CalculateFermatIterations(N);
    iterations = std::min(iterations, (size_t)(1) << 24);
    std::cout << "Trying " << iterations << " iterations of FMMod20Precomp (Fermat) factorization..." << std::endl;
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

    std::cout << "Current factors: " << factors.GetString() << std::endl;
    std::cout << "Remainder after FMMod20Precomp: " << remainder << std::endl;

    // Try using Pollards P-1
    std::cout << "Trying Pollard's P-1 (B2**32)..." << std::endl;
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

    std::cout << "Current factors: " << factors.GetString() << std::endl;
    std::cout << "Remainder after P-1: " << remainder << std::endl;

    // Use Pollard's rho algorithm
    std::cout << "Trying Brent-Pollard's rho..." << std::endl;
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

    std::cout << "Current factors: " << factors.GetString() << std::endl;
    std::cout << "Remainder after Brent-Pollard's rho: " << remainder << std::endl;

    std::cout << "Falling back to trial division..." << std::endl;
    result = TrialDivision(remainder, threads, 0, false, 0, T(modulus * threads), T(0), modulus);
    if (result) {
        factors.Update(result.value());
        return factors;
    }

    return std::nullopt;
}

}

#endif // FACTORISE_HPP