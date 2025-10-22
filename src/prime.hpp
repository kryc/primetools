#ifndef PRIME_HPP
#define PRIME_HPP

#include <vector>

#include <gmpxx.h>

#include "trial_division_data.hpp"
#include "util.hpp"

namespace primetools {

template <typename T>
const T
GetNthPrime(
    const T N
)
{
    if (N == 0) {
        throw std::invalid_argument("N must be greater than 0");
    } else if (N == 1) {
        return 2;
    } else if (N == 2) {
        return 3;
    } else if (N == 3) {
        return 5;
    }

    T candidate = 1;
    T index = 4; // Compensate for starting at 1 and skipping 2

    // Use wheel30 to skip non-prime candidates
    uint32_t gapword = primetools::WHEEL30GAPSUINT32;
    do {
        primetools::increment(candidate, gapword & 0xf);
        gapword = std::rotr(gapword, 4);

        if (primetools::isprime(candidate)) {
            index++;
        }
    } while (index <= N);

    return candidate;
}

template <typename T>
std::vector<T>
GetPrimesInRange(
    const T Lower,
    const T Upper
)
{
    if (Lower >= Upper) {
        return {};
    }

    std::vector<T> primes;

    if (Lower <= 2 && Upper >= 2) {
        primes.push_back(2);
    }
    if (Lower <= 3 && Upper >= 3) {
        primes.push_back(3);
    }
    if (Lower <= 5 && Upper >= 5) {
        primes.push_back(5);
    }

    T candidate = Lower;

    // Ensure candidate is odd
    if (candidate % 2 == 0) {
        candidate++;
    }

    // Seek backwards until candidate is congruent to 1 mod 30
    while (candidate % 30 != 1) {
        candidate--;
    }

    uint32_t gapword = primetools::WHEEL30GAPSUINT32;
    do {
        if (candidate >= Lower && primetools::isprime(candidate)) {
            primes.push_back(candidate);
        }
        primetools::increment(candidate, gapword & 0xf);
        gapword = std::rotr(gapword, 4);
    } while (candidate <= Upper);

    return primes;
}

template <typename T>
std::vector<T>
GetPrimesTo(
    const T Upper
)
{
    return GetPrimesInRange<T>(1, Upper);
}

template <typename T>
std::vector<T>
GetFirstNPrimes(
    const T N
)
{
    if (N == 0) {
        throw std::invalid_argument("N must be greater than 0");
    }

    std::vector<T> primes;
    primes.reserve(N);

    if (N >= 1) {
        primes.push_back(2);
    }
    if (N >= 2) {
        primes.push_back(3);
    }
    if (N >= 3) {
        primes.push_back(5);
    }

    if (N <= 3) {
        return primes;
    }

    T candidate = 1;

    // Use wheel30 to skip non-prime candidates
    uint32_t gapword = primetools::WHEEL30GAPSUINT32;
    do {
        primetools::increment(candidate, gapword & 0xf);
        gapword = std::rotr(gapword, 4);

        if (primetools::isprime(candidate)) {
            primes.push_back(candidate);
        }
    } while (primes.size() < N);

    return primes;
}

template <typename T>
const uint64_t
GetTrialDivisionModuliForPrime(
    const T Prime
)
{
    if (Prime < 2) {
        throw std::invalid_argument("Prime must be greater than 1");
    }

    auto primes = primetools::GetPrimesTo<T>(Prime);
    uint64_t modulus = 1;
    for (const auto& p : primes) {
        modulus *= p;
    }
    return modulus;
}

} // namespace primetools

#endif // PRIME_HPP