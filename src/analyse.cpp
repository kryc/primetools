#include <bit>

#include <gmpxx.h>

#include "analyse.hpp"

namespace primetools
{

const size_t
GuessSizeOfPrimeFactors(
    const mpz_class& N,
    const bool RoundUp
)
{
    // Calculate the size of N's constituent primes
    size_t bits = mpz_sizeinbase(N.get_mpz_t(), 2) / 2;

    // If the size is slightly smaller than a power of two, bump it up
    // ie, if bits is 63, we bump it to 64, 1023 to 1024, etc.
    if (RoundUp && !std::has_single_bit(bits)) {
        // Test if adding 1 to bits makes it a power of two
        if (std::has_single_bit(bits + 1)) {
            // If it is a power of two, just use it
            bits += 1;
        }
    }

    return bits;
}

} // namespace primetools