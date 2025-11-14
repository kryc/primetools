#ifndef BIGINT_HPP
#define BIGINT_HPP

#include <array>
#include <cstddef>
#include <algorithm>
#include <span>

#include <gmpxx.h>

template <size_t Bits>
class BigInt {
private:
    static_assert(Bits % 32 == 0, "Bits must be a multiple of 32");
    static constexpr size_t NumWords = Bits / 32;
    std::array<uint32_t, NumWords> data = {0};

    size_t limbs(void) const {
        for (size_t i = NumWords; i-- > 0;) {
            if (data[i] != 0) return i + 1;
        }
        return 0;
    }

    size_t bit_length(void) const {
        for (size_t i = NumWords; i-- > 0;) {
            if (data[i] != 0) {
                uint32_t w = data[i];
                size_t bl = 32 - __builtin_clz(w);
                return i * 32 + bl;
            }
        }
        return 0;
    }

public:
    BigInt() = default;

    // Access individual 32-bit words
    uint32_t& operator[](size_t index) {
        return data[index];
    }

    const uint32_t& operator[](size_t index) const {
        return data[index];
    }

    // Get the size in bits
    static constexpr size_t size_in_bits() {
        return Bits;
    }

    // operator =
    BigInt& operator=(const BigInt& other) {
        for (size_t i = 0; i < NumWords; ++i) {
            data[i] = other.data[i];
        }
        return *this;
    }
    BigInt& operator=(const uint64_t other) {
        data[0] = static_cast<uint32_t>(other);
        data[1] = static_cast<uint32_t>(other >> 32);
        for (size_t i = 2; i < NumWords; ++i) {
            data[i] = 0;
        }
        return *this;
    }
    BigInt& operator=(const mpz_class& other) {
        mpz_class temp = other;
        for (size_t i = 0; i < NumWords; ++i) {
            data[i] = static_cast<uint32_t>(temp.get_ui());
            temp >>= 32;
        }
        return *this;
    }

    // operator +=
    BigInt& operator+=(const uint64_t other) {
        size_t carry = other;
        for (size_t i = 0; i < NumWords; ++i) {
            uint64_t sum = static_cast<uint64_t>(data[i]) + carry;
            data[i] = static_cast<uint32_t>(sum);
            carry = sum >> 32;
            if (carry == 0) {
                break;
            }
        }
        return *this;
    }

    // operator >>=
    BigInt& operator>>=(const size_t shift) {
        if (shift >= Bits) {
            for (size_t i = 0; i < NumWords; ++i) {
                data[i] = 0;
            }
            return *this;
        }

        size_t wordShift = shift / 32;
        size_t bitShift = shift % 32;

        if (wordShift > 0) {
            for (size_t i = 0; i < NumWords - wordShift; ++i) {
                data[i] = data[i + wordShift];
            }
            for (size_t i = NumWords - wordShift; i < NumWords; ++i) {
                data[i] = 0;
            }
        }

        if (bitShift > 0) {
            for (size_t i = 0; i < NumWords - 1; ++i) {
                data[i] = (data[i] >> bitShift) | (data[i + 1] << (32 - bitShift));
            }
            data[NumWords - 1] >>= bitShift;
        }

        return *this;
    }

    // operator <<=
    BigInt& operator<<=(const size_t shift) {
        if (shift >= Bits) {
            for (size_t i = 0; i < NumWords; ++i) {
                data[i] = 0;
            }
            return *this;
        }

        size_t wordShift = shift / 32;
        size_t bitShift = shift % 32;

        if (wordShift > 0) {
            for (size_t i = NumWords; i-- > 0;) {
                if (i >= wordShift) {
                    data[i] = data[i - wordShift];
                }
            }
            for (size_t i = 0; i < wordShift; ++i) data[i] = 0;
        }

        if (bitShift > 0) {
            for (size_t i = NumWords - 1; i > 0; --i) {
                data[i] = (data[i] << bitShift) | (data[i - 1] >> (32 - bitShift));
            }
            data[0] <<= bitShift;
        }

        return *this;
    }

    // operator==
    bool operator==(const BigInt& other) const {
        for (size_t i = 0; i < NumWords; ++i) {
            if (data[i] != other.data[i]) {
                return false;
            }
        }
        return true;
    }
    bool operator==(const uint64_t other) const {
        if (data[0] != static_cast<uint32_t>(other) || data[1] != static_cast<uint32_t>(other >> 32)) {
            return false;
        }
        for (size_t i = 2; i < NumWords; ++i) {
            if (data[i] != 0) {
                return false;
            }
        }
        return true;
    }
    bool operator==(const mpz_class& other) const {
        mpz_class temp;
        for (size_t i = 0; i < NumWords; ++i) {
            temp += mpz_class(data[i]) << (32 * i);
        }
        return temp == other;
    }

    // Convert to mpz_class
    mpz_class to_mpz() const {
        mpz_class result;
        for (size_t i = 0; i < NumWords; ++i) {
            if (data[i] != 0) {
                result += mpz_class(data[i]) << (32 * i);
            }
        }
        return result;
    }

    // Check zero
    bool isZero() const {
        for (size_t i = 0; i < NumWords; ++i) {
            if (data[i] != 0) return false;
        }
        return true;
    }

    // restoring_divides
    // Check if this BigInt is evenly divisible by `divisor` using restoring division.
    bool restoring_divides(const BigInt& divisor) const {
        if (divisor.isZero()) return false;
        if (this->isZero()) return true; // 0 divisible by any non-zero divisor

        BigInt<Bits> remainder; // starts at 0
        size_t bl = bit_length();
        for (size_t bit = bl; bit-- > 0;) {
            // remainder <<= 1; remainder += current bit of dividend
            for (size_t j = NumWords - 1; j > 0; --j) {
                remainder[j] = (remainder[j] << 1) | (remainder[j - 1] >> 31);
            }
            remainder[0] = (remainder[0] << 1) | ((data[bit / 32] >> (bit % 32)) & 1u);

            // if remainder >= divisor then remainder -= divisor
            int cmp = 0;
            for (size_t j = NumWords; j-- > 0;) {
                uint32_t r = remainder[j];
                uint32_t d = divisor[j];
                if (r < d) { cmp = -1; break; }
                if (r > d) { cmp = 1; break; }
            }
            if (cmp >= 0) {
                uint64_t borrow = 0;
                for (size_t j = 0; j < NumWords; ++j) {
                    uint64_t a = remainder[j];
                    uint64_t b = divisor[j];
                    uint64_t sub = a - b - borrow;
                    remainder[j] = static_cast<uint32_t>(sub);
                    borrow = (a < b + borrow) ? 1ull : 0ull;
                }
            }
        }
        return remainder.isZero();
    }

    // divides
    // Check if this BigInt is evenly divisible by `divisor` using Knuth's Algorithm D (multi-limb division).
    bool divides(const BigInt& divisor) const {
        if (divisor.isZero()) return false;
        if (this->isZero()) return true;

        size_t n = divisor.limbs();
        size_t m_plus_n = this->limbs();
        size_t m = m_plus_n - n;

        if (n == 1) {
            uint64_t d = divisor.data[0];
            if (d == 0) return false;
            uint64_t rem = 0;
            for (size_t idx = m_plus_n; idx-- > 0;) {
                rem = (rem << 32) | data[idx];
                rem %= d;
            }
            return rem == 0;
        }

        // Copy limbs into vectors u (size m+n+1) and v (size n)
        std::array<uint32_t, Bits * 2 + 1> u = {0};
        std::array<uint32_t, Bits> v = {0};
        uint32_t* const u32 = u.data();
        uint32_t* const v32 = v.data();
        std::copy_n(data.begin(), m_plus_n, u32);
        std::copy_n(divisor.data.begin(), n, v32);

        // Normalize: left-shift so that top limb of v has high bit set
        uint32_t vn1 = v32[n - 1];
        unsigned s = 0;
        if (vn1 != 0) {
            s = static_cast<unsigned>(__builtin_clz(vn1)); // number of leading zeros in top limb
        }
        if (s > 0) {
            // shift v left by s
            uint32_t carry = 0;
            for (size_t i = 0; i < n; ++i) {
                uint64_t val = (static_cast<uint64_t>(v32[i]) << s) | carry;
                v32[i] = static_cast<uint32_t>(val);
                carry = static_cast<uint32_t>(val >> 32);
            }
            // shift u left by s
            carry = 0;
            for (size_t i = 0; i <= m + n; ++i) {
                uint64_t val = (static_cast<uint64_t>(u32[i]) << s) | carry;
                u32[i] = static_cast<uint32_t>(val);
                carry = static_cast<uint32_t>(val >> 32);
            }
            if (m + n + 1 < u.size()) {
                u32[m + n + 1] = carry;
            }
        }

        const uint64_t B = 1ull << 32;
        // Main loop j = m down to 0
        for (size_t jj = m + 1; jj-- > 0; ) {
            size_t j = jj; // current position
            // Estimate qhat
            uint64_t ujn = static_cast<uint64_t>(u32[j + n]);
            uint64_t ujn1 = static_cast<uint64_t>(u32[j + n - 1]);
            uint64_t numerator = (ujn << 32) | ujn1;
            uint64_t v1 = static_cast<uint64_t>(v32[n - 1]);
            uint64_t qhat = numerator / v1;
            uint64_t rhat = numerator % v1;
            if (qhat > 0xFFFFFFFFull) qhat = 0xFFFFFFFFull;

            if (n >= 2) {
                uint64_t v2 = static_cast<uint64_t>(v32[n - 2]);
                uint64_t ujn2 = static_cast<uint64_t>(u32[j + n - 2]);
                while (qhat * v2 > (rhat << 32) + ujn2) {
                    --qhat;
                    rhat += v1;
                    if (rhat >= B) break;
                }
            }

            // Multiply and subtract qhat * v from u segment starting at j
            uint64_t borrow = 0;
            for (size_t i = 0; i < n; ++i) {
                uint64_t p = qhat * static_cast<uint64_t>(v32[i]);
                uint64_t low = static_cast<uint32_t>(p);
                uint64_t temp = low + borrow;
                uint32_t temp_low = static_cast<uint32_t>(temp);
                uint64_t temp_carry = temp >> 32;
                uint64_t uval = static_cast<uint64_t>(u32[j + i]);
                bool sub_borrow = uval < temp_low;
                u32[j + i] = static_cast<uint32_t>(uval - temp_low);
                borrow = (p >> 32) + temp_carry + (sub_borrow ? 1ull : 0ull);
            }
            uint64_t u_top = static_cast<uint64_t>(u32[j + n]);
            bool negative = u_top < borrow;
            u32[j + n] = static_cast<uint32_t>(u_top - borrow);

            if (negative) {
                // Add back v to the segment and adjust
                uint64_t carry = 0;
                for (size_t i = 0; i < n; ++i) {
                    uint64_t sum = static_cast<uint64_t>(u32[j + i]) + static_cast<uint64_t>(v32[i]) + carry;
                    u32[j + i] = static_cast<uint32_t>(sum);
                    carry = sum >> 32;
                }
                u32[j + n] = static_cast<uint32_t>(static_cast<uint64_t>(u32[j + n]) + carry);
                // qhat would be decremented, but quotient is unused
            }

        }

        // Unnormalize remainder by right-shifting u[0..n-1] by s
        if (s > 0) {
            uint32_t prev = 0;
            for (size_t i = n; i-- > 0;) {
                uint32_t cur = u32[i];
                uint32_t newval = (cur >> s) | (prev << (32 - s));
                u32[i] = newval;
                prev = cur & ((1u << s) - 1u);
            }
        }

        // Check remainder zero
        for (size_t i = 0; i < n; ++i) {
            if (u32[i] != 0) {
                return false;
            }
        }
        return true;
    }
};

#ifdef __AVX512F__
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

    // Add a scalar uint64_t to lane 0 of all limbs with carry propagation
    BigIntAVX& operator+=(uint64_t other) {
        uint64_t carry = other;
        for (size_t w = 0; w < NumWords && carry; ++w) {
            alignas(64) uint32_t lane_vals[NumLanes];
            _mm512_store_si512(reinterpret_cast<__m512i*>(lane_vals), limbs[w]);
            uint64_t full = static_cast<uint64_t>(lane_vals[0]) + carry;
            lane_vals[0] = static_cast<uint32_t>(full);
            carry = full >> 32;
            limbs[w] = _mm512_load_si512(reinterpret_cast<const __m512i*>(lane_vals));
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

#endif // __AVX512F__

# endif // BIGINT_HPP