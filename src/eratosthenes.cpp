#include <math.h>
#include <vector>

#include "eratosthenes.hpp"
#include "euler.hpp"

std::vector<size_t>
SieveOfEratosthenes(
    const size_t Limit
)
{
    std::vector<bool> primes(Limit + 1, true);
    primes[0] = primes[1] = false; // 0 and 1 are not prime numbers

    const size_t s = static_cast<size_t>(sqrt(Limit));

    for (size_t p = 2; p < s; ++p) {
        if (primes[p]) {
            for (size_t multiple = p * p; multiple <= Limit; multiple += p) {
                primes[multiple] = false;
            }
        }
    }

    std::vector<size_t> result;
    for (size_t num = 2; num <= Limit; ++num) {
        if (primes[num]) {
            result.push_back(num);
        }
    }

    return result;
}

std::vector<size_t>
SieveOfEratosthenesQuadraticResidueP(
    const size_t Limit,
    const mpz_class P
)
{
    std::vector<bool> primes(Limit + 1, true);
    primes[0] = primes[1] = false; // 0 and 1 are not prime numbers

    const size_t s = static_cast<size_t>(sqrt(Limit));

    for (size_t p = 2; p <= s; ++p) {
        if (primes[p]) {
            // Set all multiples of p to false
            for (size_t multiple = p * p; multiple <= Limit; multiple += p) {
                primes[multiple] = false;
            }
            // Check for quadratic residue of p against P
            if (!EulerCriterion(p, P)) {
                primes[p] = false;
            }
        }
    }

    std::vector<size_t> result;
    for (size_t num = 2; num <= Limit; ++num) {
        if (primes[num]) {
            result.push_back(num);
        }
    }

    return result;
}