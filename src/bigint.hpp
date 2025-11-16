#ifndef BIGINT_HPP
#define BIGINT_HPP

#include <array>
#include <cstddef>
#include <algorithm>
#include <span>

#include <gmpxx.h>

#ifdef __AVX512F__
#include "bigint_avx.hpp"
#endif // __AVX512F__

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

# endif // BIGINT_HPP