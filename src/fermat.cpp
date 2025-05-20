#include "fermat.hpp"

namespace primetools {

const std::optional<std::pair<mpz_class, mpz_class>>
FermatFactorization(
    const mpz_class& N,
    const size_t Max
)
{
    if (N < 2) {
        return std::nullopt;
    }

    mpz_class a, b;

    a = N - 1;
    mpz_sqrt(a.get_mpz_t(), a.get_mpz_t());
    a = a + 1;

    for (size_t i = 0; i < Max; ++i) {
        b = (a * a) - N;

        if (mpz_perfect_square_p(b.get_mpz_t())) {
            return std::make_pair(a - b, a + b);
        }

        ++a;
    }

    return std::nullopt;
}

}