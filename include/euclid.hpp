#ifndef EUCLID_HPP
#define EUCLID_HPP

#include "gmpxx.h"

// Euclidean algorithm to compute the greatest
// common divisor (GCD) of two numbers
template <typename T>
static  const T
EuclideanAlgorithm(
    const T& a,
    const T& b
)
{
    T x = a;
    T y = b;

    while (y != 0) {
        const T temp = y;
        y = x % y;
        x = temp;
    }

    return x;
}

#endif // EUCLID_HPP