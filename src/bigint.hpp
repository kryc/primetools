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
    static_assert(Bits % 64 == 0, "Bits must be a multiple of 64");
    static constexpr size_t NumWords = Bits / 64;
    std::array<uint64_t, NumWords> data = {0};

    size_t limbs(void) const {
        for (size_t i = NumWords; i-- > 0;) {
            if (data[i] != 0) {
                return i + 1;
            }
        }
        return 0;
    }

    size_t bit_length(void) const {
        for (size_t i = NumWords; i-- > 0;) {
            if (data[i] != 0) {
                uint64_t w = data[i];
                size_t bl = 64 - __builtin_clzll(w);
                return i * 64 + bl;
            }
        }
        return 0;
    }

public:
    BigInt() = default;

    // Access individual 64-bit words
    uint64_t& operator[](size_t index) {
        return data[index];
    }

    const uint64_t& operator[](size_t index) const {
        return data[index];
    }

    // Get the size in bits
    static constexpr size_t size_in_bits() {
        return Bits;
    }

    // operator =
    BigInt& operator=(const BigInt& other) {
        std::copy(other.data.begin(), other.data.end(), data.begin());
        return *this;
    }
    BigInt& operator=(const uint64_t other) {
        data[0] = other;
        std::fill(data.begin() + 1, data.end(), 0);
        return *this;
    }
    BigInt& operator=(const mpz_class& other) {
        mpz_class temp = other;
        for (size_t i = 0; i < NumWords; ++i) {
            data[i] = static_cast<uint64_t>(temp.get_ui());
            temp >>= 64;
        }
        return *this;
    }

    // operator +=
    BigInt& operator+=(const uint64_t other) {
        uint64_t carry = other;
        for (size_t i = 0; i < NumWords && carry; ++i) {
            uint64_t old = data[i];
            uint64_t sum = old + carry;
            // overflow if sum < old
            carry = (sum < old) ? 1u : 0u;
            data[i] = sum;
        }
        return *this;
    }

    // operator -=
    BigInt& operator-=(const BigInt& other) {
        uint64_t borrow = 0;
        for (size_t i = 0; i < NumWords; ++i) {
            uint64_t a = data[i];
            uint64_t b = other.data[i];

            uint64_t sub = a - b;
            uint64_t borrow1 = (a < b) ? 1u : 0u;

            uint64_t sub2 = sub - borrow;
            uint64_t borrow2 = (sub < borrow) ? 1u : 0u;

            data[i] = sub2;
            borrow = borrow1 | borrow2;
        }
        return *this;
    }
    BigInt& operator-=(const uint64_t other) {
        uint64_t borrow = other;
        for (size_t i = 0; i < NumWords && borrow; ++i) {
            uint64_t a = data[i];
            uint64_t sub = a - borrow;
            uint64_t new_borrow = (a < borrow) ? 1u : 0u;
            data[i] = sub;
            borrow = new_borrow;
        }
        return *this;
    }

    // operator >>=
    BigInt& operator>>=(const size_t shift) {
        if (shift >= Bits) {
            for (size_t i = 0; i < NumWords; ++i) data[i] = 0;
            return *this;
        }

        size_t wordShift = shift / 64;
        size_t bitShift = shift % 64;

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
                data[i] = (data[i] >> bitShift) | (data[i + 1] << (64 - bitShift));
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

        size_t wordShift = shift / 64;
        size_t bitShift = shift % 64;

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
                data[i] = (data[i] << bitShift) | (data[i - 1] >> (64 - bitShift));
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
        if (data[0] != other) return false;
        for (size_t i = 1; i < NumWords; ++i) {
            if (data[i] != 0) return false;
        }
        return true;
    }
    bool operator==(const mpz_class& other) const {
        mpz_class temp;
        for (size_t i = 0; i < NumWords; ++i) {
            temp += mpz_class(data[i]) << (64 * i);
        }
        return temp == other;
    }

    // Convert to mpz_class
    mpz_class to_mpz() const {
        mpz_class result;
        for (size_t i = 0; i < NumWords; ++i) {
            if (data[i] != 0) {
                result += mpz_class(data[i]) << (64 * i);
            }
        }
        return result;
    }

    template <typename T>
    T to() const {
        if constexpr (std::is_same_v<T, mpz_class>) {
            return to_mpz();
        } else if constexpr (std::is_same_v<T, uint64_t>) {
            return data[0];
        } else if constexpr (std::is_same_v<T, __uint128_t>) {
            __uint128_t result = data[0];
            if (NumWords > 1) {
                result |= static_cast<__uint128_t>(data[1]) << 64;
            }
            return result;
        } else {
            static_assert(sizeof(T) <= 16, "Conversion to type T not supported");
        }
    }

    // Check zero
    bool zero() const {
        for (size_t i = 0; i < NumWords; ++i) {
            if (data[i] != 0) {
                return false;
            }
        }
        return true;
    }

    // restoring_divides
    // Check if this BigInt (the divisor) evenly divides `dividend` using restoring division.
    bool restoring_divides(const BigInt& dividend) const {
        // A zero divisor does not divide any value.
        if (this->zero()) return false;
        // Zero dividend is divisible by any non-zero divisor.
        if (dividend.zero()) return true;

        BigInt<Bits> remainder; // starts at 0
        for (size_t bit = dividend.bit_length(); bit-- > 0;) {
            // remainder <<= 1; remainder += current bit of dividend
            remainder <<= 1;
            remainder += (dividend.data[bit / 64] >> (bit % 64)) & 1u;

            // if remainder >= *this (divisor) then remainder -= *this
            int cmp = 0;
            for (size_t j = NumWords; j-- > 0;) {
                uint64_t r = remainder[j];
                uint64_t d = data[j];
                if (r < d) { cmp = -1; break; }
                if (r > d) { cmp = 1; break; }
            }
            if (cmp >= 0) {
                remainder -= *this;
            }
        }
        return remainder.zero();
    }

    // divides
    // Check if this BigInt (the divisor) evenly divides `dividend` using Knuth's Algorithm D (multi-limb division).
    bool divides(const BigInt& dividend) const {
        // A zero divisor does not divide any value.
        if (this->zero()) return false;
        // Zero dividend is divisible by any non-zero divisor.
        if (dividend.zero()) return true;

        size_t n = this->limbs();       // divisor length
        size_t u_len = dividend.limbs(); // dividend length
        size_t m = u_len - n;

        if (n == 1) {
            uint64_t d = data[0];
            if (d == 0) return false;

            __uint128_t rem = 0;
            for (size_t idx = u_len; idx-- > 0;) {
                rem = (rem << 64) | dividend.data[idx];
                rem %= d;
            }
            return rem == 0;
        }

        // Copy limbs into vectors u (size m+n+1) and v (size n)
        std::array<uint64_t, Bits / 64 * 2 + 2> u = {0}; // dividend
        std::array<uint64_t, Bits / 64 + 1> v = {0};     // divisor
        std::copy(dividend.data.begin(), dividend.data.begin() + u_len, u.data());
        std::copy(data.begin(), data.begin() + n, v.data());

        // Normalize: left-shift so that top limb of v has high bit set
        uint64_t vn1 = v[n - 1];
        unsigned s = 0;
        if (vn1 != 0) {
            s = static_cast<unsigned>(__builtin_clzll(vn1)); // number of leading zeros in top limb
        }
        if (s > 0) {
            // shift v left by s
            uint64_t carry = 0;
            for (size_t i = 0; i < n; ++i) {
                __uint128_t val = (static_cast<__uint128_t>(v[i]) << s) | carry;
                v[i] = static_cast<uint64_t>(val);
                carry = static_cast<uint64_t>(val >> 64);
            }
            // shift u left by s
            carry = 0;
            for (size_t i = 0; i <= m + n; ++i) {
                __uint128_t val = (static_cast<__uint128_t>(u[i]) << s) | carry;
                u[i] = static_cast<uint64_t>(val);
                carry = static_cast<uint64_t>(val >> 64);
            }
            if (m + n + 1 < u.size()) {
                u[m + n + 1] = carry;
            }
        }

        const __uint128_t B = static_cast<__uint128_t>(1) << 64;
        // Main loop j = m down to 0
        for (size_t jj = m + 1; jj-- > 0; ) {
            size_t j = jj; // current position
            // Estimate qhat
            __uint128_t ujn = static_cast<__uint128_t>(u[j + n]);
            __uint128_t ujn1 = static_cast<__uint128_t>(u[j + n - 1]);
            __uint128_t numerator = (ujn << 64) | ujn1;
            __uint128_t v1 = static_cast<__uint128_t>(v[n - 1]);
            __uint128_t qhat = numerator / v1;
            __uint128_t rhat = numerator % v1;
            if (qhat > 0xFFFFFFFFFFFFFFFFull) qhat = 0xFFFFFFFFFFFFFFFFull;

            if (n >= 2) {
                __uint128_t v2 = static_cast<__uint128_t>(v[n - 2]);
                __uint128_t ujn2 = static_cast<__uint128_t>(u[j + n - 2]);
                while (qhat * v2 > (rhat << 64) + ujn2) {
                    --qhat;
                    rhat += v1;
                    if (rhat >= B) break;
                }
            }

            // Multiply and subtract qhat * v from u segment starting at j
            __uint128_t borrow = 0;
            for (size_t i = 0; i < n; ++i) {
                __uint128_t p = qhat * static_cast<__uint128_t>(v[i]);
                __uint128_t low = static_cast<__uint128_t>(static_cast<uint64_t>(p));
                __uint128_t temp = low + borrow;
                uint64_t temp_low = static_cast<uint64_t>(temp);
                __uint128_t temp_carry = temp >> 64;
                __uint128_t uval = static_cast<__uint128_t>(u[j + i]);
                bool sub_borrow = uval < temp_low;
                u[j + i] = static_cast<uint64_t>(uval - temp_low);
                borrow = (p >> 64) + temp_carry + (sub_borrow ? 1 : 0);
            }
            __uint128_t u_top = static_cast<__uint128_t>(u[j + n]);
            bool negative = u_top < borrow;
            u[j + n] = static_cast<uint64_t>(u_top - borrow);

            if (negative) {
                // Add back v to the segment and adjust
                __uint128_t carry = 0;
                for (size_t i = 0; i < n; ++i) {
                    __uint128_t sum = static_cast<__uint128_t>(u[j + i]) + static_cast<__uint128_t>(v[i]) + carry;
                    u[j + i] = static_cast<uint64_t>(sum);
                    carry = sum >> 64;
                }
                u[j + n] = static_cast<uint64_t>(static_cast<__uint128_t>(u[j + n]) + carry);
                // qhat would be decremented, but quotient is unused
            }

        }

        // Unnormalize remainder by right-shifting u[0..n-1] by s
        if (s > 0) {
            uint64_t prev = 0;
            for (size_t i = n; i-- > 0;) {
                uint64_t cur = u[i];
                uint64_t newval = (cur >> s) | (prev << (64 - s));
                u[i] = newval;
                prev = cur & ((uint64_t(1) << s) - 1u);
            }
        }

        // Check remainder zero
        for (size_t i = 0; i < n; ++i) {
            if (u[i] != 0) {
                return false;
            }
        }
        return true;
    }
};

# endif // BIGINT_HPP