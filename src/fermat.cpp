#include "fermat.hpp"
#include <algorithm>
#include <iostream>


namespace primetools {

const std::optional<std::pair<mpz_class, mpz_class>>
FermatFactorisation(
    const mpz_class& N,
    const size_t Offset,
    const size_t Max
)
{
    if (N < 2) {
        return std::nullopt;
    }

    mpz_class a, b;

    mpz_sqrt(a.get_mpz_t(), N.get_mpz_t());
    a += 1 + Offset;

    for (size_t i = 0; i < Max; ++i) {
        b = (a * a) - N;

        if (mpz_perfect_square_p(b.get_mpz_t())) {
            mpz_sqrt(b.get_mpz_t(), b.get_mpz_t());
            return std::make_pair(a - b, a + b);
        }

        ++a;
    }

    return std::nullopt;
}


// Sieve improvement for fermat's factorization
const std::optional<std::pair<mpz_class, mpz_class>>
FermatSieve(
    const mpz_class& N,
    const size_t Offset,
    const size_t Max
)
{
    if (N < 2) {
        return std::nullopt;
    }

    mpz_class a, b;

    mpz_sqrt(a.get_mpz_t(), N.get_mpz_t());
    a += 1 + Offset;

    for (size_t i = 0; i < Max; ++i) {

        // Sieve improvement: Skip values that cannot be perfect squares
        // We do it as a switch statement to help the compiler optimize
        auto remainder = mpz_fdiv_ui(a.get_mpz_t(), 20);

        bool skip = false;
        switch (remainder)
        {
        case 1:
        case 9:
        case 11:
        case 19:
            skip = true;
            break;
        }

        if (!skip){
            b = (a * a) - N;
            if (mpz_perfect_square_p(b.get_mpz_t())) {
                mpz_sqrt(b.get_mpz_t(), b.get_mpz_t());
                return std::make_pair(a - b, a + b);
            }
        }

        ++a;
    }

    return std::nullopt;
}

}