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

    // Access specific lane and limb (word index)
    uint32_t lane_word(size_t lane, size_t word) const {
        alignas(64) uint32_t tmp[NumLanes];
        _mm512_store_si512(reinterpret_cast<__m512i*>(tmp), limbs[word]);
        return tmp[lane];
    }

    // Assign each lane from span<uint64_t> as a little-endian integer.
    // src[lane] is interpreted as a 64-bit value for that lane; remaining
    // higher limbs stay zero.
    BigIntAVX& operator=(std::span<const uint64_t> src) {
        // zero all limbs/lanes
        for (size_t w = 0; w < NumWords; ++w) {
            limbs[w] = _mm512_setzero_si512();
        }
        if (src.empty()) return *this;
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

    // Assign from mpz_class in span
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

    // Assign a single mpz_class to all lanes
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

    // Add a scalar uint32_t to all lanes with carry propagation across 32-bit limbs.
    BigIntAVX& operator+=(uint32_t other) {
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

    bool restoring_divides(const BigIntAVX& divisor) const {
        std::array<__m512i, NumWords> remainder{};
        const __m512i zero = _mm512_setzero_si512();
        const __m512i one  = _mm512_set1_epi32(1);

        // Fast path: if any lane of the divisor is zero, treat as not divisible
        // (mirrors scalar BigInt::restoring_divides behaviour on zero divisor).
        __mmask16 any_divisor_nonzero = 0;
        for (size_t j = 0; j < NumWords; ++j) {
            any_divisor_nonzero |= _mm512_cmp_epi32_mask(divisor.limbs[j], zero, _MM_CMPINT_NE);
        }
        // If any lane has divisor == 0, bail out.
        if (any_divisor_nonzero != 0xFFFF) {
            return false;
        }

        for (size_t bit = max_bitlength(); bit-- > 0;) {
            // Shift remainder left by 1 across all lanes
            for (size_t j = NumWords - 1; j > 0; --j) {
                remainder[j] = _mm512_or_si512(_mm512_slli_epi32(remainder[j], 1),
                                              _mm512_srli_epi32(remainder[j - 1], 31));

            }
            remainder[0] = _mm512_or_si512(_mm512_slli_epi32(remainder[0], 1),
                                          _mm512_and_si512(_mm512_srli_epi32(limbs[0], bit % 32),
                                                           _mm512_set1_epi32(1)));

            // if remainder >= divisor then remainder -= divisor
            // Robust multi-limb comparison per lane using an "undecided" mask.
            __mmask16 lt_mask_all = 0;      // lanes where remainder < divisor
            __mmask16 gt_mask_all = 0;      // lanes where remainder > divisor (kept for clarity/debug, not used)
            __mmask16 undecided   = 0xFFFF; // lanes still equal so far
            for (size_t j = NumWords; j-- > 0;) {
                __m512i r = remainder[j];
                __m512i d = divisor.limbs[j];
                __mmask16 lt_mask = _mm512_cmp_epi32_mask(r, d, _MM_CMPINT_LT) & undecided;
                __mmask16 gt_mask = _mm512_cmp_epi32_mask(r, d, _MM_CMPINT_GT) & undecided;

                lt_mask_all |= lt_mask;
                gt_mask_all |= gt_mask;
                undecided   &= ~(lt_mask | gt_mask);
            }
            // lanes with remainder >= divisor: not (remainder < divisor)
            (void)gt_mask_all; // suppress unused warning
            __mmask16 ge_mask = ~lt_mask_all;

            // Subtract divisor from remainder on those lanes, propagating borrow across limbs.
            __mmask16 borrow_mask = 0;
            for (size_t j = 0; j < NumWords; ++j) {
                __m512i r = remainder[j];
                __m512i d = divisor.limbs[j];

                // r - d
                __m512i diff = _mm512_sub_epi32(r, d);
                // Detect per-lane borrow: r < d
                __mmask16 this_borrow = _mm512_cmp_epi32_mask(r, d, _MM_CMPINT_LT);

                // Apply existing borrow into r: if borrow_mask set, subtract 1
                __m512i borrow_vec = _mm512_mask_blend_epi32(borrow_mask, zero, one);
                diff = _mm512_sub_epi32(diff, borrow_vec);

                // New borrow is: (r < d) OR (r == d AND borrow_in)
                __mmask16 eq_mask = _mm512_cmp_epi32_mask(r, d, _MM_CMPINT_EQ);
                __mmask16 new_borrow = this_borrow | (eq_mask & borrow_mask);

                // Keep diff only on lanes where we actually subtract (ge_mask)
                remainder[j] = _mm512_mask_mov_epi32(r, ge_mask, diff);

                borrow_mask = new_borrow & ge_mask;
            }
        }
        // Compare all lanes with zero across all limbs.
        uint16_t zero_mask = 0xFFFF;
        for (size_t j = 0; j < NumWords; ++j) {
            __mmask16 cmp = _mm512_cmpeq_epi32_mask(remainder[j], zero);
            // Keep only lanes that are still zero across all limbs.
            zero_mask &= static_cast<uint16_t>(cmp);
        }
        // Any lane divisible iff at least one lane's remainder is zero.
        return zero_mask != 0x0000;
    }
};

#endif // BIGINT_AVX_HPP