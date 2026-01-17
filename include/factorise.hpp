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
    size_t modulus = 510510;
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
    if (primetools::isprime(remainder)) {
        factors.AddFactor(remainder);
        return factors;
    } else if (remainder == 1) {
        return factors;
    }

    Log("Current factors: " + factors.GetString(), Verbose);
    Log("Remainder after FMMod20Precomp: " + remainder.get_str(), Verbose);

    // Try using Pollards P-1
    Log("Trying Pollard's P-1...", Verbose);
    for (;;) {
        resultpair = PollardsPMinus1MT(remainder, threads);
        if (resultpair) {
            // If either factor is prime, add it to the factors
            if (primetools::isprime(resultpair->first)) {
                factors.AddFactor(resultpair->first);
            }
            if (primetools::isprime(resultpair->second)) {
                factors.AddFactor(resultpair->second);
            }

            remainder = N / factors.Product();
            if (primetools::isprime(remainder)) {
                factors.AddFactor(remainder);
                return factors;
            } else if (remainder == 1) {
                return factors;
            }
            Log("Current factors: " + factors.GetString(), Verbose);
        } else {
            break;
        }
    }

    remainder = N / factors.Product();
    if (primetools::isprime(remainder)) {
        factors.AddFactor(remainder);
        return factors;
    } else if (remainder == 1) {
        return factors;
    }

    Log("Current factors: " + factors.GetString(), Verbose);
    Log("Remainder after P-1: " + remainder.get_str(), Verbose);

    // Use Pollard's rho algorithm
    Log("Trying Brent-Pollard's rho...", Verbose);
    for (;;) {
        resultpair = BrentPollardsRhoMT(remainder, threads, DefaultM, (size_t)(1) << 4);
        if (resultpair) {
            // If either factor is prime, add it to the factors
            if (primetools::isprime(resultpair->first)) {
                factors.AddFactor(resultpair->first);
            }
            if (primetools::isprime(resultpair->second)) {
                factors.AddFactor(resultpair->second);
            }

            remainder = N / factors.Product();
            if (primetools::isprime(remainder)) {
                factors.AddFactor(remainder);
                return factors;
            } else if (remainder == 1) {
                return factors;
            }
            Log("Current factors: " + factors.GetString(), Verbose);
        } else {
            break;
        }        
    }

    // Fall back to trial division if still not fully factored
    remainder = N / factors.Product();
    if (primetools::isprime(remainder)) {
        factors.AddFactor(remainder);
        return factors;
    } else if (remainder == 1) {
        return factors;
    }

    Log("Current factors: " + factors.GetString(), Verbose);
    Log("Remainder after Brent-Pollard's rho: " + remainder.get_str(), Verbose);

    Log("Falling back to trial division...", Verbose);
    T last_finish = (modulus * threads);
    modulus = 9699690;
    result = TrialDivision(remainder, threads, 0, false, 0, last_finish, T(0), modulus);
    if (result) {
        factors.Update(result.value());
        return factors;
    }

    return std::nullopt;
}

}

#endif // FACTORISE_HPP