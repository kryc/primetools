#ifndef EULER_HPP
#define EULER_HPP

#include <gmpxx.h>

const mpz_class
EulerTotient(
    const mpz_class N
);

const bool
EulerCriterion(
    const mpz_class A,
    const mpz_class P
);

#endif // EULER_HPP