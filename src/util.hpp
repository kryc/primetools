#ifndef Umpz_classIL_HPP
#define Umpz_classIL_HPP

#include <bit>
#include <cmath>
#include <stdexcept>
#include <stdint.h>
#include <gmp.h>
#include <gmpxx.h>

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

const mpz_class
floor_div(
    const mpz_class& A,
    const mpz_class& B
);

const bool
very_fast_prime_test(
    const mpz_class& Candidate
);

const bool
fast_prime_test(
    const mpz_class& Candidate
);

static inline
const bool
divides(
    const mpz_class& N,
    const mpz_class& Value
)
{
    return mpz_divisible_p(N.get_mpz_t(), Value.get_mpz_t());
}

static inline
const bool
divides(
    const mpz_class& N,
    const uint64_t Value
)
{
    return mpz_divisible_ui_p(N.get_mpz_t(), Value);
}

static inline
const bool
divides(
    const uint64_t N,
    const uint64_t Value
)
{
    return (N % Value) == 0;
}

static inline void
increment(
    mpz_class& Value,
    const mpz_class& Amount = 1
)
{
    mpz_add(Value.get_mpz_t(), Value.get_mpz_t(), Amount.get_mpz_t());
}

static inline void
increment(
    mpz_class& Value,
    const uint64_t Amount = 1
)
{
    mpz_add_ui(Value.get_mpz_t(), Value.get_mpz_t(), Amount);
}

static inline void
increment(
    mpz_class& Value,
    const uint32_t Amount = 1
)
{
    mpz_add_ui(Value.get_mpz_t(), Value.get_mpz_t(), Amount);
}

// Cover all basic types
template <typename T1, typename T2>
static inline void
increment(
    T1& Value,
    const T2 Amount = 1
)
{
    Value += Amount;
}

static inline void
toggle_bit(
    mpz_class& Value,
    const size_t Bit
)
{
    mpz_combit(Value.get_mpz_t(), Bit);
}

template <typename T>
static inline void
toggle_bit(
    T& Value,
    const size_t Bit
)
{
    Value ^= (T(1) << Bit);
}

static inline
const size_t
bit_size(
    const mpz_class& Value
)
{
    return mpz_sizeinbase(Value.get_mpz_t(), 2);
}

template <typename T>
static inline
const size_t
bit_size(
    const T& Value
)
{
    return std::bit_width(Value);
}

} // namespace primetools

#endif