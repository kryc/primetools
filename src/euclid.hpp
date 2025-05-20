#ifndef EUCLID_HPP
#define EUCLID_HPP

#include "gmpxx.h"

// Euclidean algorithm to compute the greatest
// common divisor (GCD) of two numbers
inline const mpz_class
EuclideanAlgorithm(
    const mpz_class& a,
    const mpz_class& b
)
{
    mpz_class x = a;
    mpz_class y = b;

    while (y != 0) {
        mpz_class temp = y;
        y = x % y;
        x = temp;
    }

    return x;
}

#endif // EUCLID_HPP