#include <bit>
#include <cassert>
#include <iostream>
#include <optional>

#ifdef __AVX512F__
#include <immintrin.h>
#endif

#include <gmpxx.h>

#include "factors.hpp"
#include "trial_division.hpp"

namespace primetools {

const std::optional<PrimeFactors<mpz_class>>
TrialDivisionSimd(
    const mpz_class& N,
    const size_t MaxIterations
)
{
#ifdef __AVX512F__
    if (N < 2) {
        return std::nullopt;
    }

    if (!N.fits_ulong_p() || N.get_ui() > (1UL << 53)) {
        std::cout << "Number is too large for SIMD trial division." << std::endl;
        return std::nullopt;
    }

    uint64_t target64 = N.get_ui();

    // Get the square root of N as this is our upper bound for trial division
    const mpz_class upper = sqrt(N);
    // const double upperd = upper.get_d();

    std::cout << "Trying factorization of primes below " << upper.get_str() << " using SIMD wheel-30 factorization." << std::endl;

    // Initialize the first 16 candidates
    const __m512i targetu64 = _mm512_set1_epi64(target64);
    // Convert to a double
    const __m512d target = _mm512_cvtepu64_pd(targetu64);
    // Initialize the candidates and step
    __m512i candidates = _mm512_set_epi64(7, 11, 13, 17, 19, 23, 29, 31);
    const __m512i step = _mm512_set1_epi64(30);
    const __m512i zero = _mm512_set1_epi64(0);

    for (size_t i = 0; i < MaxIterations; i += 1)
    {
        // Compute the approximate quotients
        const __m512d candidatesd = _mm512_cvtepu64_pd(candidates);
        const __m512d quotients = _mm512_div_pd(target, candidatesd);

        // Round the quotients to the nearest integer
        const __m512d rounded = _mm512_roundscale_pd(quotients, _MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC);

        // Convert the quotients back to integers
        const __m512i rounded64 = _mm512_cvttpd_epu64(rounded);

        // Multiply back to get the remainders
        const __m512i products = _mm512_mullox_epi64(rounded64, candidates);

        // Calculate the remainders
        const __m512i remainders = _mm512_sub_epi64(targetu64, products);

        // Check if any of the remainders are zero
        const __mmask8 mask = _mm512_cmpeq_epi64_mask(remainders, zero);
        
        if (mask != 0) {
            int index = __builtin_ctzll(mask);

            std::cout << "Found a factor!" << std::endl;
            alignas(64) uint64_t factor_array[8];
            _mm512_store_epi64(factor_array, candidates);
            const uint64_t factor = factor_array[index];
            if (mpz_divisible_ui_p(N.get_mpz_t(), factor)) {
                return std::make_pair(mpz_class(factor), N / factor);
            }
        }
        // Increment the candidates by the step
        candidates = _mm512_add_epi64(candidates, step);
    }
#else
    std::cerr << "SIMD trial division is not supported on this platform." << std::endl;
#endif
    return std::nullopt;
}

} // namespace primetools