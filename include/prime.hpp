#ifndef PRIME_HPP
#define PRIME_HPP

#include <array>
#include <span>
#include <vector>

#include <gmpxx.h>

#include "maths.hpp"
#include "primegenerator.hpp"
#include "util.hpp"
#include "wheel30.hpp"

namespace primetools {

const std::span<const uint64_t>
GetSmallPrimes(
    void
);

template <typename T>
const T
GetNthPrime(
    const T N
)
{
    if (N == 0) {
        throw std::invalid_argument("N must be greater than 0");
    }
    if (N <= 100) {
        return GetSmallPrimes()[N - 1];
    }

    T candidate = 1;
    T index = 4; // Compensate for starting at 1 and skipping 2

    // Use wheel30 to skip non-prime candidates
    uint32_t gapword = kWheel30;
    do {
        primetools::increment((uint64_t&)candidate, gapword & kWheel30Mask);
        gapword = std::rotr(gapword, kWheel30BitsPerGap);

        if (primetools::IsPrime(candidate)) {
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
    PrimeGenerator<T> generator(9699690);

    for (;;) {
        const T prime = generator.Next();
        if (prime >= Lower && prime <= Upper) {
            primes.push_back(prime);
        }
        if (prime >= Upper) {
            break;
        }
    }

    return primes;
}

std::span<const uint64_t>
GetPrimesTo(
    const uint64_t Upper
);

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
    PrimeGenerator<T> generator;
    primes.reserve(N);

    for (T i = 0; i < N; ++i) {
        primes.push_back(generator.Next());
    }

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

    auto primes = primetools::GetPrimesTo(Prime);
    uint64_t modulus = 1;
    for (const auto& p : primes) {
        modulus *= p;
    }
    return modulus;
}

const bool LoadPrimeGaps(
    std::string_view FilePath
);

const bool LoadPrimes(
    std::string_view FilePath
);

void LoadPrimeGapsInNewThread(
    std::string_view FilePath
);

const size_t
GetCachedPrimesLimit(
    void
);

} // namespace primetools

#endif // PRIME_HPP