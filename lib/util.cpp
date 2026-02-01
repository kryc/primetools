#include <cctype>
#include <iomanip>
#include <sstream>
#include <string>


#include "gmp.h"
#include "gmpxx.h"

#include "miller_rabin.hpp"
#include "util.hpp"

namespace primetools {

const bool IsPrime(
    const mpz_class& N
)
{
    // Select a K for the Miller-Rabin test based on N
    return mpz_probab_prime_p(N.get_mpz_t(), 20) > 0;
}

const bool IsPrime(
    const __uint128_t N
)
{
    return MillerRabin(N, 20);
}

const bool IsPrime(
    const uint64_t N
)
{
    return MillerRabin(N, 10);
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