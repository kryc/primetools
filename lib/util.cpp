#include <cctype>
#include <iomanip>
#include <sstream>
#include <string>


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
    return mpz_probab_prime_p(N.get_mpz_t(), 20) > 0;
}

const bool isprime(
    const uint64_t N
)
{
    mpz_class n(static_cast<unsigned long>(N));
    return isprime(n);
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

const mpz_class
floor_div(
    const mpz_class& A,
    const mpz_class& B)
{
    mpz_class res;
    mpz_fdiv_q(res.get_mpz_t(), A.get_mpz_t(), B.get_mpz_t());
    return res;
}

const bool
very_fast_prime_test(
    const mpz_class& Candidate
)
{
    const mpz_srcptr candidate = Candidate.get_mpz_t();

    if (
        mpz_divisible_ui_p(candidate, 3) ||
        mpz_divisible_ui_p(candidate, 5) ||
        mpz_divisible_ui_p(candidate, 7)
    )
    {
        return false;
    }

    return true;
}

const bool
fast_prime_test(
    const mpz_class& Candidate
)
{
    const mpz_srcptr candidate = Candidate.get_mpz_t();

    if (very_fast_prime_test(Candidate) == false) {
        return false; // Quick check for small primes
    }

    // Trial division with small primes
    int small_primes[] = {11, 13, 17, 19, 23, 29};
    for (int p : small_primes) {
        if (mpz_divisible_ui_p(candidate, p))
        {
            return false;
        }
    }

    // Perform a reduced Miller-Rabin test (only 5 rounds)
    return mpz_probab_prime_p(candidate, 5) > 0;
}

const bool
is_numeric(
    const std::string_view Str
)
{
    if (Str.empty()) {
        return false;
    }

    for (size_t i = 0; i < Str.size(); ++i) {
        if (!std::isdigit(Str[i])) {
            return false;
        }
    }

    return true;
}

const std::string
GetHex(
    const uint32_t Value,
    const size_t Width
)
{
    std::stringstream ss;
    ss << std::hex << std::setfill('0') << std::setw(Width) << Value;
    return ss.str();
}

}