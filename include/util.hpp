#ifndef Umpz_classIL_HPP
#define Umpz_classIL_HPP

#include <bit>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <string_view>
#include <ostream>
#include <type_traits>
#include <vector>

#include <stdint.h>
#include <gmp.h>
#include <gmpxx.h>

#include "euclid.hpp"

// Convert __uint128_t to string (base 10)
inline std::string uint128_to_string(__uint128_t value) {
    if (value == 0) return "0";

    std::string result;
    while (value > 0) {
        unsigned digit = static_cast<unsigned>(value % 10);
        result.insert(result.begin(), static_cast<char>('0' + digit));
        value /= 10;
    }
    return result;
}

// Overload operator<< for __uint128_t
inline std::ostream& operator<<(std::ostream& os, __uint128_t value) {
    std::string s = uint128_to_string(value);
    os.write(s.data(), static_cast<std::streamsize>(s.size()));
    return os;
}

// Optional: Overload for __int128_t (signed)
inline std::ostream& operator<<(std::ostream& os, __int128_t value) {
    if (value < 0) {
        os.put('-');
        value = -value;
    }
    return os << static_cast<__uint128_t>(value);
}

namespace primetools {

static inline std::string
ToString(
    const mpz_class& Number
)
{
    return Number.get_str();
}

static inline std::string
ToString(
    const __uint128_t& Number
)
{
    return uint128_to_string(Number);
}

static inline std::string
ToString(
    const uint64_t& Number
)
{
    return std::to_string(Number);
}

static inline std::string
ToSuperScriptDigit(
    const size_t Digit
)
{
    static const std::array<std::string, 10> superscript_digits = {"⁰","¹","²","³","⁴","⁵","⁶","⁷","⁸","⁹"};
    if (Digit > 9) {
        throw std::invalid_argument("Digit must be in range 0-9");
    }
    return superscript_digits[Digit];
}

template <typename T>
static inline std::string
ToSuperScript(
    const T& Number
)
{
    T n = Number;
    std::string result;
    size_t digit;
    while (n > 9) {
        if constexpr (std::is_same_v<T, mpz_class>) {
            mpz_class remainder;
            mpz_fdiv_r_ui(remainder.get_mpz_t(), n.get_mpz_t(), 10);
            digit = remainder.get_ui();
        } else {
            digit = n % 10;
        }
        result = ToSuperScriptDigit(digit) + result;
        n /= 10;
    }
    if constexpr (std::is_same_v<T, mpz_class>) {
        digit = n.get_ui();
    } else {
        digit = n;
    }
    result = ToSuperScriptDigit(digit) + result;
    return result;
}

template <typename T>
const std::string
TruncateNumber(
    const T& Number,
    const size_t StartDigits = 5,
    const size_t EndDigits = 3
)
{
    std::string str;
    if constexpr (std::is_same_v<T, __uint128_t>) {
        str = uint128_to_string(Number);
    }
    else {
        mpz_class n = Number;
        str = n.get_str();
    }
    if (str.length() <= StartDigits + EndDigits) {
        return str;
    }
    return str.substr(0, StartDigits) + "..." + str.substr(str.length() - EndDigits);
}

const bool IsPrime(
    const mpz_class& N
);

const bool IsPrime(
    const __uint128_t N
);

const bool IsPrime(
    const uint64_t N
);

const bool
very_fast_prime_test(
    const mpz_class& Candidate
);

const bool
fast_prime_test(
    const mpz_class& Candidate
);

static inline void
toggle_bit(
    mpz_class& Value,
    const size_t Bit
)
{
    mpz_combit(Value.get_mpz_t(), Bit);
}

template <typename T>
static inline void
toggle_bit(
    T& Value,
    const size_t Bit
)
{
    Value ^= (T(1) << Bit);
}

static inline
const size_t
BitSize(
    const mpz_class& Value
)
{
    return mpz_sizeinbase(Value.get_mpz_t(), 2);
}

template <typename T>
static inline
const size_t
BitSize(
    const T& Value
)
{
    return std::bit_width(Value);
}

static inline
const bool
fits_uint128(
    const mpz_class& Value
)
{
    return BitSize(Value) <= 128;
}

static inline __uint128_t
MpzToUint128(
    const mpz_class& Value
)
{
    __uint128_t result = Value.get_ui();
    __uint128_t high = mpz_class(Value >> 64).get_ui();
    result |= (high << 64);
    return result;
}

static inline mpz_class
Uint128ToMpz(
    const __uint128_t& Value
)
{
    mpz_class result = mpz_class(static_cast<uint64_t>(Value & 0xFFFFFFFFFFFFFFFF));
    result |= (mpz_class(static_cast<uint64_t>(Value >> 64)) << 64);
    return result;
}

template <typename T1, typename T2>
static inline T1
ConvertType(
    const T2& Value
)
{
    if constexpr (std::is_same_v<T1, mpz_class>) {
        if constexpr (std::is_same_v<T2, __uint128_t>) {
            return Uint128ToMpz(Value);
        } else {
            return mpz_class(Value);
        }
    } else if constexpr (std::is_same_v<T1, __uint128_t>) {
        if constexpr (std::is_same_v<T2, mpz_class>) {
            if (mpz_sgn(Value.get_mpz_t()) < 0) {
                throw std::overflow_error("ConvertType: Value negative for target type");
            }
            if (!fits_uint128(Value)) {
                throw std::overflow_error("ConvertType: Value too large for target type");
            }
            return MpzToUint128(Value);
        } else {
            return static_cast<__uint128_t>(Value);
        }
    } else {
        if (Value > std::numeric_limits<T1>::max()) {
            throw std::overflow_error("ConvertType: Value too large for target type");
        }
        if constexpr (std::is_same_v<T2, mpz_class>) {
            return static_cast<T1>(Value.get_ui());
        } else {
            return static_cast<T1>(Value);
        }
    }
}

const bool
is_numeric(
    const std::string_view Str
);

static inline
const mpz_class
min(
    const mpz_class& A,
    const mpz_class& B
)
{
    return (A < B) ? A : B;
}

static inline
const uint64_t
min(
    const uint64_t A,
    const uint64_t B
)
{
    return (A < B) ? A : B;
}

static inline
const mpz_class
max(
    const mpz_class& A,
    const mpz_class& B
)
{
    return (A > B) ? A : B;
}

static inline
const uint64_t
max(
    const uint64_t A,
    const uint64_t B
)
{
    return (A > B) ? A : B;
}

static inline
const uint64_t
get_ui(
    const mpz_class& Value
)
{
    return Value.get_ui();
}

static inline
const uint64_t
get_ui(
    const __uint128_t& Value
)
{
    return static_cast<uint64_t>(Value);
}

const std::string
GetHex(
    const uint32_t Value,
    const size_t Width = 0
);

template<typename T>
std::vector<uint8_t>
SerializeVLE(
    const T Value
)
{
    std::vector<uint8_t> bytes;
    T val = Value;
    do {
        uint8_t byte;
        if constexpr (std::is_same_v<T, mpz_class>) {
            byte = static_cast<uint8_t>(val.get_ui() & 0x7F);
        } else {
            byte = static_cast<uint8_t>(val & 0x7F);
        }
        val >>= 7;
        if (val != 0) {
            byte |= 0x80; // Set continuation bit
        }
        bytes.push_back(byte);
    } while (val != 0);
    return bytes;
}

} // namespace primetools

#endif