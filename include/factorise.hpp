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

constexpr size_t kSmallPminus1B2 = 20;
constexpr size_t kSmallPminus1B = (size_t)(1) << kSmallPminus1B2;
constexpr size_t kSmallPminus1Bases = 128;
constexpr size_t kLargePminus1B2 = 28;
constexpr size_t kLargePminus1B = (size_t)(1) << kLargePminus1B2;
constexpr size_t kLargePminus1Bases = 32;
constexpr size_t kSmallWheelModulus = 510510;
constexpr size_t kLargeWheelModulus = 6469693230;
constexpr size_t kFermatMaxIterations2 = 14;
constexpr size_t kFermatMaxIterations = (size_t)(1) << kFermatMaxIterations2;

// Factorise a perfect square
template <typename T>
const std::optional<PrimeFactors<T>>
FactorisePerfectSquare(
    const T& N
)
{
    if (N < 2) {
        return std::nullopt;
    }

    if (IsPerfectSquare(N)) {
        const T sqrtN = primetools::Sqrt(N);
        return PrimeFactors<T>::FromPair(sqrtN, sqrtN);
    }

    return std::nullopt;
}

// Validate our factorization result and remainder
template <typename T>
static inline const bool
ValidateFactorization(
    const T& N,
    const PrimeFactors<T>& Factors
)
{
    const T product = Factors.Product();
    if (product > N) {
        return false;
    }
    if (N % product != 0) {
        return false;
    }
    return true;
}

template <typename T>
static inline void
ValidateFactorizationOrThrow(
    const T& N,
    const PrimeFactors<T>& Factors
)
{
    if (!ValidateFactorization(N, Factors)) {
        throw std::runtime_error("Invalid factorization result. N = " + ToString(N) + ", Factors = " + Factors.GetString() + ", Product = " + ToString(Factors.Product()));
    }
}

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
    T remainder = N;
    size_t threads = Threads > 0 ? Threads : std::thread::hardware_concurrency();

    // Check for perfect square
    LogFn("Checking for perfect square...");
    auto result = FactorisePerfectSquare(N);
    if (result) {
        return result;
    }

    // Check for small prime factors
    LogFn("Checking small primes...");
    result = TrialDivision(N, 0, kSmallWheelModulus, false, 0, T(0), T(kSmallWheelModulus * threads), kSmallWheelModulus, TrialDivisionStrategy::Linear, LogFn);
    if (result) {
        factors = result.value();
        remainder = N / result->Product();
        if (primetools::IsPrime(remainder)) {
            factors.AddFactor(remainder);
            return factors;
        } else if (remainder == 1) {
            return factors;
        }
    }

    // Validate current factorization
    ValidateFactorizationOrThrow(N, factors);

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
    LogFn(factors.GetString() + " R: " + ToString(remainder) + " A: P-1 (B=2^" + std::to_string(kSmallPminus1B2) + ")");
    constexpr size_t kPollardB = kSmallPminus1B;
    constexpr size_t kPollardBases = kSmallPminus1Bases;
    auto resultpair = PollardsPMinus1MT(remainder, threads, kPollardB, kPollardBases);
    if (resultpair) {
        // If either factor is prime, add it to the factors
        if (primetools::IsPrime(resultpair->first)) {
            factors.AddFactor(resultpair->first);
        }
        if (primetools::IsPrime(resultpair->second)) {
            factors.AddFactor(resultpair->second);
        }
    }

    remainder = N / factors.Product();
    if (primetools::IsPrime(remainder)) {
        factors.AddFactor(remainder);
        return factors;
    } else if (remainder == 1) {
        return factors;
    }

    // Validate current factorization
    ValidateFactorizationOrThrow(N, factors);

    // Use Fermat's factorization method up to 2^24 iterations
    LogFn(factors.GetString() + " R: " + ToString(remainder) + " A: FMMod20Precomp");
    size_t iterations = CalculateFermatIterations(N);
    iterations = std::min(iterations, kFermatMaxIterations);
    resultpair = FMMod20PrecompMT(remainder, threads, iterations);
    if (resultpair) {
        // If either factor is prime, add it to the factors
        if (primetools::IsPrime(resultpair->first)) {
            factors.AddFactor(resultpair->first);
        }
        if (primetools::IsPrime(resultpair->second)) {
            factors.AddFactor(resultpair->second);
        }
    }

    remainder = N / factors.Product();
    if (primetools::IsPrime(remainder)) {
        factors.AddFactor(remainder);
        return factors;
    } else if (remainder == 1) {
        return factors;
    }

    // Validate current factorization
    ValidateFactorizationOrThrow(N, factors);

    // Use Pollard's rho algorithm
    LogFn(factors.GetString() + " R: " + ToString(remainder) + " A: Brent-Pollard's Rho");
    for (;;) {
        resultpair = BrentPollardsRhoMT(remainder, threads, DefaultM, (size_t)(1) << 4);
        if (resultpair) {
            // If either factor is prime, add it to the factors
            if (primetools::IsPrime(resultpair->first)) {
                factors.AddFactor(resultpair->first);
            }
            if (primetools::IsPrime(resultpair->second)) {
                factors.AddFactor(resultpair->second);
            }

            remainder = N / factors.Product();
            if (primetools::IsPrime(remainder)) {
                factors.AddFactor(remainder);
                return factors;
            } else if (remainder == 1) {
                return factors;
            }
            LogFn(factors.GetString());
        } else {
            break;
        }        
    }

    // // Try using Pollards P-1 with a large B
    LogFn(factors.GetString() + " R: " + ToString(remainder) + " A: P-1 (B=2^" + std::to_string(kLargePminus1B2) + ")");
    constexpr size_t kPollardLargeB = kLargePminus1B;
    constexpr size_t kPollardLargeBases = kLargePminus1Bases;
    resultpair = PollardsPMinus1MT(remainder, threads, kPollardLargeB, kPollardLargeBases);
    if (resultpair) {
        // If either factor is prime, add it to the factors
        if (primetools::IsPrime(resultpair->first)) {
            factors.AddFactor(resultpair->first);
        }
        if (primetools::IsPrime(resultpair->second)) {
            factors.AddFactor(resultpair->second);
        }
    }

    remainder = N / factors.Product();
    if (primetools::IsPrime(remainder)) {
        factors.AddFactor(remainder);
        return factors;
    } else if (remainder == 1) {
        return factors;
    }

    // Validate current factorization
    ValidateFactorizationOrThrow(N, factors);

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
    if (primetools::IsPrime(remainder)) {
        factors.AddFactor(remainder);
        return factors;
    } else if (remainder == 1) {
        return factors;
    }

    // Validate current factorization
    ValidateFactorizationOrThrow(N, factors);

    // Fall back to hybrid trial division (this is not ideal!)
    LogFn(factors.GetString() + " R: " + ToString(remainder) + " A: Trial Division");
    result = TrialDivision(remainder, threads, 0, false, 0, T(0), T(0), kLargeWheelModulus, TrialDivisionStrategy::MeetInTheMiddle, LogFn);
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
    // See if we can convert the number to a smaller type for faster factorization
    if constexpr (std::is_same_v<T, mpz_class>) {
        if (N.fits_ulong_p()) {
            auto factors = Factorise<uint64_t>((uint64_t)N.get_ui(), Threads, FactorDBPath, LogFn);
            if (factors) {
                return factors->template Convert<mpz_class>();
            } else {
                return std::nullopt;
            }
        } else if (mpz_sizeinbase(N.get_mpz_t(), 2) <= 128) {
            auto factors = Factorise<__uint128_t>(MpzToUint128(N), Threads, FactorDBPath, LogFn);
            if (factors) {
                return factors->template Convert<mpz_class>();
            } else {
                return std::nullopt;
            }
        }
    }

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

inline static void
PreCacheWheelAndPrimes(
    void
)
{
    GetWheelGapsForModulus(kLargeWheelModulus);
    GetPrimesTo(kLargePminus1B);
}

inline static void
PreCacheWheelAndPrimesInNewThread(
    void
)
{
    std::thread([]() {
        GetWheelGapsForModulus(kLargeWheelModulus);
    }).detach();
    std::thread([]() {
        GetPrimesTo(kLargePminus1B);
    }).detach();
}

} // namespace primetools

#endif // FACTORISE_HPP