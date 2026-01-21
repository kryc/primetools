#ifndef FACTORISE_HPP
#define FACTORISE_HPP

#include <algorithm>
#include <cstdint>
#include <functional>
#include <iostream>
#include <string_view>
#include <thread>
#include <optional>
#include <utility>
#include <vector>

#include <gmpxx.h>

#include "factordb.hpp"
#include "factors.hpp"
#include "fermat.hpp"
#include "logging.hpp"
#include "pollard.hpp"
#include "shanks.hpp"
#include "trial_division.hpp"

namespace primetools {

// Factorise a perfect square
const std::optional<PrimeFactors<mpz_class>>
FactorisePerfectSquare(
    const mpz_class& N
);

// Try to factorise with increasing complexity
template <typename T>
const std::optional<PrimeFactors<T>>
FactoriseNumber(
    const T& N,
    FactorDB<T>& Database,
    const size_t Threads = 0,
    LogCallback LogFn = LogStdOut
)
{
    if (N < 2) {
        return std::nullopt;
    }

    PrimeFactors<T> factors;
    size_t threads = Threads > 0 ? Threads : std::thread::hardware_concurrency();

    // Check for perfect square
    LogFn("Checking for perfect square...");
    auto result = FactorisePerfectSquare(N);
    if (result) {
        return result;
    }

    // Check for small prime factors
    LogFn("Checking small primes...");
    const size_t modulus = 510510;
    result = TrialDivision(N, 0, modulus, false, 0, T(0), T(modulus * threads), modulus, TrialDivisionStrategy::Linear, LogFn);
    if (result) {
        factors = result.value();
    }

    T remainder = N / result->Product();
    if (primetools::isprime(remainder)) {
        factors.AddFactor(remainder);
        return factors;
    } else if (remainder == 1) {
        return factors;
    }

    // Do a quick lookup of the remainder in the factor DB
    if (Database.IsOpen()) {
        auto lookup = Database.GetFactors(remainder);
        if (lookup) {
            LogFn("Found cached factors in FactorDB.");
            factors.Update(lookup.value());
            return factors;
        }
    }

    // Try using Pollards P-1 with a small B
    LogFn("F: " + factors.GetString() + " R: " + remainder.get_str() + " A: P-1 (B=2^20)");
    constexpr size_t kPollardB = (size_t)(1) << 20;
    constexpr size_t kPollardBases = 128;
    auto resultpair = PollardsPMinus1MT(remainder, threads, kPollardB, kPollardBases);
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

    // Use Fermat's factorization method up to 2^24 iterations
    LogFn("F: " + factors.GetString() + " R: " + remainder.get_str() + " A: FMMod20Precomp");
    size_t iterations = CalculateFermatIterations(N);
    iterations = std::min(iterations, (size_t)8192*2);
    resultpair = FMMod20PrecompMT(remainder, threads, iterations);
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

    // Use Pollard's rho algorithm
    LogFn("F: " + factors.GetString() + " R: " + remainder.get_str() + " A: Brent-Pollard's Rho");
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
            LogFn("F: " + factors.GetString());
        } else {
            break;
        }        
    }

    // // Try using Pollards P-1 with a large B
    // LogFn("F: " + factors.GetString() + " R: " + remainder.get_str() + " A: P-1 (B=2^30)");
    // constexpr size_t kPollardLargeB = (size_t)(1) << 30;
    // constexpr size_t kPollardLargeBases = 128;
    // resultpair = PollardsPMinus1MT(remainder, threads, kPollardLargeB, kPollardLargeBases);
    // if (resultpair) {
    //     // If either factor is prime, add it to the factors
    //     if (primetools::isprime(resultpair->first)) {
    //         factors.AddFactor(resultpair->first);
    //     }
    //     if (primetools::isprime(resultpair->second)) {
    //         factors.AddFactor(resultpair->second);
    //     }
    // }

    // remainder = N / factors.Product();
    // if (primetools::isprime(remainder)) {
    //     factors.AddFactor(remainder);
    //     return factors;
    // } else if (remainder == 1) {
    //     return factors;
    // }

    // Check the database for any cached factors before trial division
    if (Database.IsOpen()) {
        LogFn("Checking FactorDB for cached factors...");
        auto lookup = Database.GetFactors(remainder);
        if (lookup) {
            LogFn("Found cached factors in FactorDB.");
            factors.Update(lookup.value());
            return factors;
        }
    }

    remainder = N / factors.Product();
    if (primetools::isprime(remainder)) {
        factors.AddFactor(remainder);
        return factors;
    } else if (remainder == 1) {
        return factors;
    }

    // Fall back to hybrid trial division (this is not ideal!)
    LogFn("F: " + factors.GetString() + " R: " + remainder.get_str() + " A: Trial Division");
    T last_finish = (modulus * threads);
    result = TrialDivision(remainder, threads, modulus, false, 0, last_finish, T(0), modulus, TrialDivisionStrategy::Hybrid, LogFn);
    if (result) {
        factors.Update(result.value());
        return factors;
    }

    return std::nullopt;
}

// Try to factorise with increasing complexity
template <typename T>
const std::optional<PrimeFactors<T>>
Factorise(
    const T& N,
    const size_t Threads = 0,
    const std::string_view FactorDBPath = "",
    LogCallback LogFn = LogStdOut
)
{
    FactorDB<T> db(FactorDBPath);
    if (db.IsOpen()) {
        LogFn("Checking FactorDB for cached factors...");
        auto lookup = db.GetFactors(N);
        if (lookup) {
            LogFn("Found cached factors in FactorDB.");
            return lookup;
        }
    }
    auto factors = FactoriseNumber<T>(N, db, Threads, LogFn);
    if (factors && db.IsOpen()) {
        LogFn("Storing factors in FactorDB...");
        db.AddFactors(factors.value());
    }
    return factors;
}

}

#endif // FACTORISE_HPP