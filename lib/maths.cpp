#include "maths.hpp"

#include <cmath>
#include <limits>

namespace primetools {

double
ln_mpz(
    const mpz_class& N
)
{
    // Avoid infinities/NaNs: some build flags treat them as UB.
    // Callers in this codebase only use ln_mpz for N > 1.
    if (N <= 0) {
        return 0.0;
    }
    long exp2 = 0;
    const double mantissa = mpz_get_d_2exp(&exp2, N.get_mpz_t());
    // N = mantissa * 2^exp2 with mantissa in [0.5, 1).
    return std::log(mantissa) + static_cast<double>(exp2) * std::log(2.0);
}

int
legendre_symbol(
    uint64_t a,
    uint64_t p
)
{
    if (p == 2) {
        return (a & 1u) ? 1 : 0;
    }
    a %= p;
    if (a == 0) {
        return 0;
    }
    const uint64_t ls = static_cast<uint64_t>(primetools::ModExp(static_cast<__uint128_t>(a), (p - 1) / 2, p));
    if (ls == 1) {
        return 1;
    }
    if (ls == p - 1) {
        return -1;
    }
    return 0;
}

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
ModExp(
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
ModExp(
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

const mpz_class
floor_div(
    const mpz_class& A,
    const mpz_class& B)
{
    mpz_class res;
    mpz_fdiv_q(res.get_mpz_t(), A.get_mpz_t(), B.get_mpz_t());
    return res;
}

} // namespace primetools