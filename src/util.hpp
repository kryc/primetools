#ifndef Umpz_classIL_HPP
#define Umpz_classIL_HPP

#include <stdint.h>
#include <gmp.h>
#include <gmpxx.h>
#include <cmath>

#include <stdexcept>

#include "euclid.hpp"
#include "miller_rabin.hpp"

namespace primetools {

const mpz_class
divmod(
    const mpz_class& Value,
    const mpz_class& Divisor,
    mpz_class& Remainder
);

const mpz_class
modexp(
    const mpz_class& Base,
    const mpz_class& Exponent,
    const mpz_class& Modulus
);

const mpz_class
modexp(
    const mpz_class& Base,
    const int64_t Exponent,
    const mpz_class& Modulus
);

const mpz_class
abs(
    const mpz_class& Value
);

const bool isprime(
    const mpz_class& N
);

const mpz_class
gcd(
    const mpz_class& A,
    const mpz_class& B
);

} // namespace primetools

#endif