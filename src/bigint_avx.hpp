#ifndef BIGINT_AVX_HPP
#define BIGINT_AVX_HPP

#include <array>
#include <cstddef>
#include <algorithm>
#include <span>

#include <immintrin.h>

// BigIntAVX represents a batch of 16 independent big integers of size Bits,
// in a structure-of-arrays layout: for each limb index i we store a __m512i
// whose lanes are that limb across 16 integers.
template <size_t Bits>
class BigIntAVX {
private:
    static_assert(Bits % 32 == 0, "Bits must be a multiple of 32");
    static constexpr size_t NumWords = Bits / 32;
    static constexpr size_t NumLanes = 512 / 32; // 16 32-bit lanes per __m512i
    std::array<__m512i, NumWords> limbs{}; // limbs[word] has lanes for 16 integers
    
    size_t limbs_used(size_t lane) const {
        for (size_t i = NumWords; i-- > 0;) {
            alignas(64) uint32_t lane_vals[NumLanes];
            _mm512_store_si512(reinterpret_cast<__m512i*>(lane_vals), limbs[i]);
            if (lane_vals[lane] != 0) return i + 1;
        }
        return 0;
    }

public:
    BigIntAVX() {
        for (size_t i = 0; i < NumWords; ++i) {
            limbs[i] = _mm512_setzero_si512();
        }
    }

    mpz_class to_mpz(size_t lane) const {
        mpz_class result;
        for (size_t w = 0; w < NumWords; ++w) {
            alignas(64) uint32_t lane_vals[NumLanes];
            _mm512_store_si512(reinterpret_cast<__m512i*>(lane_vals), limbs[w]);
            if (lane_vals[lane] != 0) {
                result += mpz_class(lane_vals[lane]) << (32 * w);
            }
        }
        return result;
    }

    template <typename T>
    T to(const size_t lane) const {
        if constexpr (std::is_same_v<T, mpz_class>) {
            return to_mpz(lane);
        } else if constexpr (std::is_same_v<T, uint64_t>) {
            uint64_t result = 0;
            for (size_t w = 0; w < NumWords && w < 2; ++w) {
                alignas(64) uint32_t lane_vals[NumLanes];
                _mm512_store_si512(reinterpret_cast<__m512i*>(lane_vals), limbs[w]);
                result |= static_cast<uint64_t>(lane_vals[lane]) << (32 * w);
            }
            return result;
        } else if constexpr (std::is_same_v<T, __uint128_t>) {
            __uint128_t result = 0;
            for (size_t w = 0; w < NumWords && w < 4; ++w) {
                alignas(64) uint32_t lane_vals[NumLanes];
                _mm512_store_si512(reinterpret_cast<__m512i*>(lane_vals), limbs[w]);
                result |= static_cast<__uint128_t>(lane_vals[lane]) << (32 * w);
            }
            return result;
        } else {
            static_assert(sizeof(T) <= 16, "Conversion to type T not supported");
        }
    }

    // Access specific lane and limb (word index)
    uint32_t lane_word(size_t lane, size_t word) const {
        alignas(64) uint32_t tmp[NumLanes];
        _mm512_store_si512(reinterpret_cast<__m512i*>(tmp), limbs[word]);
        return tmp[lane];
    }

    BigIntAVX& operator=(std::span<const uint64_t> src) {
        std::fill(limbs.begin(), limbs.end(), _mm512_setzero_si512());
        size_t lanes_to_fill = std::min(src.size(), static_cast<size_t>(NumLanes));

        // low 32 bits limb 0
        {
            alignas(64) uint32_t lane_vals[NumLanes] = {0};
            for (size_t lane = 0; lane < lanes_to_fill; ++lane) {
                lane_vals[lane] = static_cast<uint32_t>(src[lane]);
            }
            limbs[0] = _mm512_load_si512(reinterpret_cast<const __m512i*>(lane_vals));
        }

        // high 32 bits limb 1, if present
        if (NumWords > 1) {
            alignas(64) uint32_t lane_vals[NumLanes] = {0};
            for (size_t lane = 0; lane < lanes_to_fill; ++lane) {
                lane_vals[lane] = static_cast<uint32_t>(src[lane] >> 32);
            }
            limbs[1] = _mm512_load_si512(reinterpret_cast<const __m512i*>(lane_vals));
        }
        return *this;
    }
    BigIntAVX& operator=(std::span<const mpz_class> vals) {
        for (size_t w = 0; w < NumWords; ++w) {
            limbs[w] = _mm512_setzero_si512();
        }
        size_t lanes_to_fill = std::min(vals.size(), static_cast<size_t>(NumLanes));
        for (size_t lane = 0; lane < lanes_to_fill; ++lane) {
            mpz_class temp = vals[lane];
            for (size_t w = 0; w < NumWords; ++w) {
                uint32_t limb = static_cast<uint32_t>(temp.get_ui());
                alignas(64) uint32_t lane_vals[NumLanes];
                _mm512_store_si512(reinterpret_cast<__m512i*>(lane_vals), limbs[w]);
                lane_vals[lane] = limb;
                limbs[w] = _mm512_load_si512(reinterpret_cast<const __m512i*>(lane_vals));
                temp >>= 32;
            }
        }
        return *this;
    }
    BigIntAVX& operator=(const std::span<const __uint128_t> vals) {
        for (size_t w = 0; w < NumWords; ++w) {
            limbs[w] = _mm512_setzero_si512();
        }
        size_t lanes_to_fill = std::min(vals.size(), static_cast<size_t>(NumLanes));
        for (size_t lane = 0; lane < lanes_to_fill; ++lane) {
            __uint128_t temp = vals[lane];
            for (size_t w = 0; w < NumWords; ++w) {
                uint32_t limb = static_cast<uint32_t>(temp & 0xFFFFFFFFull);
                alignas(64) uint32_t lane_vals[NumLanes];
                _mm512_store_si512(reinterpret_cast<__m512i*>(lane_vals), limbs[w]);
                lane_vals[lane] = limb;
                limbs[w] = _mm512_load_si512(reinterpret_cast<const __m512i*>(lane_vals));
                temp >>= 32;
            }
        }
        return *this;
    }
    BigIntAVX& operator=(const mpz_class& val) {
        for (size_t w = 0; w < NumWords; ++w) {
            alignas(64) uint32_t lane_vals[NumLanes];
            uint32_t limb = static_cast<uint32_t>(val.get_ui() >> (32 * w));
            for (size_t lane = 0; lane < NumLanes; ++lane) {
                lane_vals[lane] = limb;
            }
            limbs[w] = _mm512_load_si512(reinterpret_cast<const __m512i*>(lane_vals));
        }
        return *this;
    }
    BigIntAVX& operator=(const uint64_t val) {
        std::fill(limbs.begin(), limbs.end(), _mm512_setzero_si512());
        __m512i lane_vals = _mm512_set1_epi32(static_cast<int32_t>(val));
        limbs[0] = lane_vals;
        if (NumWords > 1) {
            __m512i lane_vals_high = _mm512_set1_epi32(static_cast<int32_t>(val >> 32));
            limbs[1] = lane_vals_high;
        }
        return *this;
    }

    BigIntAVX& operator+=(const uint32_t other) {
        const __m512i zero  = _mm512_setzero_si512();
        __m512i carry       = _mm512_set1_epi32(static_cast<int32_t>(other)); // initial delta=other

        for (size_t w = 0; w < NumWords; ++w) {
            if (w != 0 && _mm512_test_epi32_mask(carry, carry) == 0) {
                break; // no carry left to propagate
            }

            const __m512i limb = limbs[w];
            const __m512i sum  = _mm512_add_epi32(limb, carry);
            const __mmask16 carry_mask = _mm512_cmp_epu32_mask(sum, limb, _MM_CMPINT_LT);

            limbs[w] = sum;
            carry = _mm512_mask_set1_epi32(zero, carry_mask, 1); // next delta = carry (0/1)
        }
        return *this;
    }

    BigIntAVX& operator-=(const uint32_t other) {
        const __m512i zero  = _mm512_setzero_si512();
        __m512i borrow      = _mm512_set1_epi32(static_cast<int32_t>(other)); // initial delta=other

        for (size_t w = 0; w < NumWords; ++w) {
            if (w != 0 && _mm512_test_epi32_mask(borrow, borrow) == 0) {
                break; // no borrow left to propagate
            }

            const __m512i limb = limbs[w];
            const __m512i diff  = _mm512_sub_epi32(limb, borrow);
            const __mmask16 borrow_mask = _mm512_cmp_epu32_mask(limb, diff, _MM_CMPINT_LT);

            limbs[w] = diff;
            borrow = _mm512_mask_set1_epi32(zero, borrow_mask, 1);
        }
        return *this;
    }
    BigIntAVX& operator-=(const BigIntAVX& other) {
        const __m512i zero  = _mm512_setzero_si512();
        __m512i borrow      = _mm512_setzero_si512(); // initial delta=0

        for (size_t w = 0; w < NumWords; ++w) {
            const __m512i a = limbs[w];
            const __m512i b = other.limbs[w];

            // a - b
            __m512i diff = _mm512_sub_epi32(a, b);
            // Detect per-lane borrow: a < b
            __mmask16 this_borrow = _mm512_cmp_epu32_mask(a, b, _MM_CMPINT_LT);

            // Apply existing borrow into a: if borrow_mask set, subtract 1
            __m512i borrow_vec = _mm512_mask_set1_epi32(zero, this_borrow, 1);
            diff = _mm512_sub_epi32(diff, borrow_vec);

            // New borrow is: (a < b) OR (a == b AND previous borrow)
            __mmask16 equal_mask = _mm512_cmp_epu32_mask(a, b, _MM_CMPINT_EQ);
            __mmask16 new_borrow_mask = this_borrow | (equal_mask & _mm512_cmp_epu32_mask(borrow, zero, _MM_CMPINT_NE));

            limbs[w] = diff;
            borrow = _mm512_mask_set1_epi32(zero, new_borrow_mask, 1); // next delta = borrow (0/1)
        }
        return *this;
    }

    BigIntAVX& operator<<=(size_t shift) {
        if (shift == 0) return *this;
        size_t limb_shift = shift / 32;
        size_t bit_shift = shift % 32;

        if (limb_shift >= NumWords) {
            // Shift exceeds size, zero all
            for (size_t w = 0; w < NumWords; ++w) {
                limbs[w] = _mm512_setzero_si512();
            }
            return *this;
        }

        if (bit_shift == 0) {
            // Limb-aligned shift
            for (size_t w = NumWords; w-- > limb_shift;) {
                limbs[w] = limbs[w - limb_shift];
            }
        } else {
            // Bit shift with carry between limbs
            for (size_t w = NumWords; w-- > limb_shift + 1;) {
                __m512i high_part = _mm512_srli_epi32(limbs[w - limb_shift - 1], 32 - bit_shift);
                __m512i low_part  = _mm512_slli_epi32(limbs[w - limb_shift], bit_shift);
                limbs[w] = _mm512_or_si512(low_part, high_part);
            }
            // Handle the first shifted limb separately
            limbs[limb_shift] = _mm512_slli_epi32(limbs[0], bit_shift);
        }

        // Zero out the lower limbs
        for (size_t w = 0; w < limb_shift; ++w) {
            limbs[w] = _mm512_setzero_si512();
        }

        return *this;
    }

    size_t max_bitlength() const {
        // Find the highest bit set in any lane using AVX512-only logic.
        const __m512i zero = _mm512_setzero_si512();
        const __m512i top_bit_mask = _mm512_set1_epi32(1u << 31);

        // Scan limbs from most-significant to least-significant.
        for (size_t limb = NumWords; limb-- > 0;) {
            __m512i val = limbs[limb];
            // Skip limb entirely if all lanes are zero.
            if (_mm512_cmp_epi32_mask(val, zero, _MM_CMPINT_EQ) == 0xFFFF) {
                continue;
            }

            // Scan bits within this limb, from most-significant to least.
            for (int bit = 31; bit >= 0; --bit) {
                __mmask16 cmp = _mm512_cmp_epi32_mask(_mm512_and_si512(val, top_bit_mask), zero, _MM_CMPINT_NE);
                if (cmp != 0) {
                    // Global bit index is limb*32 + (bit position within limb) + 1.
                    size_t bit_in_limb = static_cast<size_t>(bit);
                    return limb * 32 + bit_in_limb + 1;
                }
                // Shift left so next bit moves into the top position.
                val = _mm512_slli_epi32(val, 1);
            }
        }
        return 0;
    }

    template <typename T>
    const bool restoring_divides(const T& Dividend) const {
        constexpr size_t DividendMaxBits  = 8192;
        constexpr size_t DividendMaxWords = DividendMaxBits / 32;

        BigIntAVX<Bits> remainder; // starts at 0
        const __m512i zero = _mm512_setzero_si512();

        // Break the dividend into an array of simd u32 limbs (LSW first).
        size_t dividend_bits = primetools::bit_size(Dividend);
        size_t dividend_words = (dividend_bits + 31) / 32;

        std::array<uint32_t, DividendMaxWords> dividend_limbs{};
        std::array<__m512i, DividendMaxWords> dividend_limbs_vec{};
        
        if constexpr (std::is_same_v<T, mpz_class>) {
            size_t count = 0;
            mpz_export(dividend_limbs.data(), &count, -1, sizeof(uint32_t), 0, 0, Dividend.get_mpz_t());
            if (count != dividend_words) {
                throw std::runtime_error("mpz_export count does not match expected dividend words");
            }
        } else {
            for (size_t w = 0; w < dividend_words; ++w) {
                dividend_limbs[w] = static_cast<uint32_t>((Dividend >> (32 * w)) & 0xFFFFFFFFull);
            }
        }
        
        // Load dividend limbs into __m512i vectors
        for (size_t w = 0; w < dividend_words; ++w) {
            dividend_limbs_vec[w] = _mm512_set1_epi32(dividend_limbs[w]);
        }

        // Reject the case where all lanes of the divisor are zero.
        __mmask16 non_zero_mask = 0;
        for (size_t j = 0; j < NumWords; ++j) {
            non_zero_mask |= _mm512_cmp_epi32_mask(limbs[j], zero, _MM_CMPINT_NE);
        }
        if (non_zero_mask == 0) {
            return false;
        }

        // Restoring division: does this divisor divide dividend in any lane?
        for (size_t bit = dividend_bits; bit-- > 0;) {
            // remainder <<= 1
            remainder <<= 1;

            // remainder += current bit of dividend (same bit for all lanes),
            // with proper multi-limb carry propagation.
            const size_t   word_idx    = bit / 32;
            const uint32_t bit_in_word = static_cast<uint32_t>(bit % 32);

            __m512i word_vec  = dividend_limbs_vec[word_idx];
            __m512i shifted   = _mm512_srli_epi32(word_vec, bit_in_word);
            __m512i ones      = _mm512_set1_epi32(1);
            __m512i masked    = _mm512_and_si512(shifted, ones);
            __mmask16 add_mask = _mm512_cmp_epu32_mask(masked, zero, _MM_CMPINT_NE);

            // Propagate a +1 across limbs where the bit is set.
            __m512i carry_vec = _mm512_mask_set1_epi32(zero, add_mask, 1);
            for (size_t j = 0; j < NumWords; ++j) {
                if (j > 0) {
                    // Stop if there is no carry left in any lane.
                    if (_mm512_test_epi32_mask(carry_vec, carry_vec) == 0) break;
                }
                __m512i limb = remainder.limbs[j];
                __m512i sum  = _mm512_add_epi32(limb, carry_vec);
                // Unsigned carry detection: sum < limb
                __mmask16 new_carry_mask = _mm512_cmp_epu32_mask(sum, limb, _MM_CMPINT_LT);
                remainder.limbs[j] = sum;
                carry_vec = _mm512_mask_set1_epi32(zero, new_carry_mask, 1);
            }

            // if remainder >= divisor then remainder -= divisor (per lane)
            __mmask16 lt_mask_all = 0;
            __mmask16 gt_mask_all = 0;
            __mmask16 undecided   = non_zero_mask;
            for (size_t j = NumWords; j-- > 0;) {
                const __m512i r = remainder.limbs[j];
                const __m512i d = limbs[j];
                const __mmask16 lt_mask = _mm512_cmp_epu32_mask(r, d, _MM_CMPINT_LT) & undecided;
                const __mmask16 gt_mask = _mm512_cmp_epu32_mask(r, d, _MM_CMPINT_GT) & undecided;
                lt_mask_all |= lt_mask;
                gt_mask_all |= gt_mask;
                undecided   &= ~(lt_mask | gt_mask);
            }
            (void)gt_mask_all; // suppress unused warning

            __mmask16 ge_mask = non_zero_mask & ~lt_mask_all; // lanes where remainder >= divisor

            // Subtract divisor from remainder on those lanes, propagating borrow across limbs.
            __mmask16 borrow_mask = 0;
            for (size_t j = 0; j < NumWords; ++j) {
                const __m512i r = remainder.limbs[j];
                const __m512i d = limbs[j];
                __m512i diff    = _mm512_sub_epi32(r, d);

                // Detect per-lane borrow: r < d (unsigned)
                const __mmask16 this_borrow = _mm512_cmp_epu32_mask(r, d, _MM_CMPINT_LT);

                // Apply existing borrow into r: if borrow_mask set, subtract 1
                const __m512i borrow_vec = _mm512_mask_set1_epi32(zero, borrow_mask, 1);
                diff = _mm512_sub_epi32(diff, borrow_vec);

                // New borrow is: (r < d) OR (r == d AND borrow_in)
                const __mmask16 eq_mask    = _mm512_cmp_epu32_mask(r, d, _MM_CMPINT_EQ);
                const __mmask16 new_borrow = this_borrow | (eq_mask & borrow_mask);

                // Keep diff only on lanes where we actually subtract (ge_mask)
                remainder.limbs[j] = _mm512_mask_mov_epi32(r, ge_mask, diff);

                borrow_mask = new_borrow & ge_mask;
            }
        }

        // Compare all lanes with zero across all limbs.
        uint16_t zero_mask = 0xFFFF;
        for (size_t j = 0; j < NumWords; ++j) {
            __mmask16 cmp = _mm512_cmpeq_epi32_mask(remainder.limbs[j], zero);
            zero_mask &= static_cast<uint16_t>(cmp);
        }

        // We consider division successful if any lane with a non-zero divisor has zero remainder.
        return (zero_mask & non_zero_mask) != 0;
    }

    template <typename T>
    const bool divides(const T& Dividend) const {
        return restoring_divides(Dividend);
    }
};

#endif // BIGINT_AVX_HPP