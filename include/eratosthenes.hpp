#ifndef ERATOSTHENES_HPP
#define ERATOSTHENES_HPP

#include <vector>

#include <gmpxx.h>

std::vector<size_t>
SieveOfEratosthenes(
    const size_t Limit
);

std::vector<size_t>
SieveOfEratosthenesQuadraticResidueP(
    const size_t Limit,
    const mpz_class P
);

#endif // ERATOSTHENES_HPP