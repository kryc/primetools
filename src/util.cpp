#include "gmp.h"
#include "gmpxx.h"

#include "util.hpp"

namespace primetools {

const mpz_class
divmod(
    const mpz_class& Value,
    const mpz_class& Divisor,
    mpz_class& Remainder
)
{
    mpz_class quotient;
    mpz_fdiv_qr(
        quotient.get_mpz_t(),
        Remainder.get_mpz_t(),
        Value.get_mpz_t(),
        Divisor.get_mpz_t()
    );
    return quotient;
}

const mpz_class
modexp(
    const mpz_class& Base,
    const mpz_class& Exponent,
    const mpz_class& Modulus
)
{
    mpz_class result;
    mpz_powm(
        result.get_mpz_t(),
        Base.get_mpz_t(),
        Exponent.get_mpz_t(),
        Modulus.get_mpz_t()
    );
    return result;
}

const mpz_class
modexp(
    const mpz_class& Base,
    const int64_t Exponent,
    const mpz_class& Modulus
)
{
    mpz_class result;
    mpz_powm_ui(
        result.get_mpz_t(),
        Base.get_mpz_t(),
        Exponent,
        Modulus.get_mpz_t()
    );
    return result;
}

const mpz_class
abs(
    const mpz_class& Value
)
{
    if (mpz_sgn(Value.get_mpz_t()) < 0) {
        mpz_class result;
        mpz_abs(result.get_mpz_t(), Value.get_mpz_t());
        return result;
    }
    return Value;
}

const bool isprime(
    const mpz_class& N
)
{
    // Select a K for the Miller-Rabin test based on N
    return MillerRabin(N, 40);
}

const mpz_class
gcd(
    const mpz_class& A,
    const mpz_class& B
)
{
    mpz_class result;
    mpz_gcd(result.get_mpz_t(), A.get_mpz_t(), B.get_mpz_t());
    return result;
}

}