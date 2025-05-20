#ifndef MILLERRABIN_HPP
#define MILLERRABIN_HPP

#include <gmpxx.h>

const bool
MillerRabin(
    const mpz_class N,
    const size_t K
);

#endif // MILLERRABIN_HPP