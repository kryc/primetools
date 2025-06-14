#ifndef ANALYSE_HPP
#define ANALYSE_HPP

#include <cstddef>

#include <gmpxx.h>

namespace primetools
{

const size_t
GuessSizeOfPrimeFactors(
    const mpz_class& N,
    const bool RoundUp = true
);

}

#endif