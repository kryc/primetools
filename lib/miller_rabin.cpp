#include "gmpxx.h"
#include <stdint.h>

#include "miller_rabin.hpp"
#include "util.hpp"

const bool
MillerRabin(
    const mpz_class N,
    const size_t K
)
{
    if (N <= 1 || N == 4) return false;
    if (N <= 3) return true;

    // Find d such that N = 2^r * d + 1 with d odd
    mpz_class d = N - 1;
    int r = 0;
    while (d % 2 == 0) {
        d /= 2;
        r++;
    }

    // Witness loop
    for (int i = 0; i < K; i++) {
        mpz_class a = rand() % (N - 4) + 2; // Random base in [2, N-2]
        mpz_class x = primetools::modexp(a, d, N);
        if (x == 1 || x == N - 1) continue;

        bool composite = true;
        for (int j = 0; j < r - 1; j++) {
            x = primetools::modexp(x, 2, N);
            if (x == N - 1) {
                composite = false;
                break;
            }
        }
        if (composite) return false;
    }
    return true;
}