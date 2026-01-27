#ifndef ERATOSTHENES_HPP
#define ERATOSTHENES_HPP

#include <cmath>
#include <vector>

#include <gmpxx.h>

namespace primetools {

static const std::vector<bool>
SieveOfEratosthenesBool(
    const size_t Limit
)
{
    std::vector<bool> primes(Limit + 1, true);
    primes[0] = primes[1] = false; // 0 and 1 are not prime numbers

    const size_t s = static_cast<size_t>(std::sqrt(Limit));

    for (size_t p = 2; p < s; ++p) {
        if (primes[p]) {
            for (size_t multiple = p * p; multiple <= Limit; multiple += p) {
                primes[multiple] = false;
            }
        }
    }

    return primes;
}

static const std::vector<bool>
SieveOfEratosthenesOddBool(
    const size_t Limit
)
{
    std::vector<bool> primes((Limit / 2) + 1, true); // Only store odd numbers
    primes[0] = false; // 1 is not a prime number

    const size_t s = static_cast<size_t>(std::sqrt(Limit));

    for (size_t p = 3; p <= s; p += 2) {
        if (primes[p / 2]) {
            for (size_t multiple = p * p; multiple <= Limit; multiple += 2 * p) {
                primes[multiple / 2] = false;
            }
        }
    }

    return primes;
}

template <typename T>
static const std::vector<T>
SieveOfEratosthenes(
    const size_t Limit
)
{
    std::vector<bool> is_prime = SieveOfEratosthenesOddBool(Limit);
    std::vector<T> primes = {2};

    for (size_t index = 0; index < is_prime.size(); ++index) {
        if (is_prime[index]) {
            primes.push_back(static_cast<T>(index * 2 + 1));
        }
    }

    return primes;
}

template <typename T>
static const std::vector<T>
SegmentedSieve(
    const size_t Limit
)
{
    const size_t segment_size = std::max(static_cast<size_t>(std::sqrt(Limit)), static_cast<size_t>(32768));
    std::vector<T> primes;
    primes.reserve(Limit / std::log(Limit)); // Approximate number of primes

    std::vector<bool> is_prime(segment_size);

    for (size_t low = 2; low <= Limit; low += segment_size) {
        std::fill(is_prime.begin(), is_prime.end(), true);
        size_t high = std::min(low + segment_size - 1, Limit);

        for (size_t p = 2; p * p <= high; ++p) {
            if (std::find(primes.begin(), primes.end(), static_cast<T>(p)) != primes.end() || p == 2) {
                size_t start = std::max(p * p, (low + p - 1) / p * p);
                for (size_t multiple = start; multiple <= high; multiple += p) {
                    is_prime[multiple - low] = false;
                }
            }
        }

        for (size_t i = low; i <= high; ++i) {
            if (is_prime[i - low]) {
                primes.push_back(static_cast<T>(i));
            }
        }
    }

    return primes;
}

template <typename T, typename T2>
static const std::vector<T>
SieveOfEratosthenesQuadraticResidueP(
    const size_t Limit,
    const T2& P
)
{
    std::vector<bool> is_prime = SieveOfEratosthenesOddBool(Limit);
    std::vector<T> primes = {2};

    for (size_t index = 0; index < is_prime.size(); ++index) {
        if (is_prime[index]) {
            T candidate = static_cast<T>(index * 2 + 1);
            if (!EulerCriterion(candidate, P)) {
                primes.push_back(candidate);
            }
        }
    }

    return primes;
}

} // namespace primetools

#endif // ERATOSTHENES_HPP