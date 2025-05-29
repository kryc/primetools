#include <math.h>
#include <vector>

#include "eratosthenes.hpp"

std::vector<unsigned long int>
SieveOfEratosthenes(
    unsigned long int limit
)
{
    std::vector<bool> is_prime(limit + 1, true);
    is_prime[0] = is_prime[1] = false; // 0 and 1 are not prime numbers

    unsigned long int sqrt_limit = static_cast<unsigned long int>(sqrt(limit));

    for (unsigned long int p = 2; p < sqrt_limit; ++p) {
        if (is_prime[p]) {
            for (unsigned long int multiple = p * p; multiple <= limit; multiple += p) {
                is_prime[multiple] = false;
            }
        }
    }

    std::vector<unsigned long int> primes;
    for (unsigned long int num = 2; num <= limit; ++num) {
        if (is_prime[num]) {
            primes.push_back(num);
        }
    }

    return primes;
}