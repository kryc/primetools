#ifndef Umpz_classIL_HPP
#define Umpz_classIL_HPP

#include <bit>
#include <cmath>
#include <stdexcept>
#include <string>
#include <string_view>
#include <ostream>
#include <vector>

#include <stdint.h>
#include <gmp.h>
#include <gmpxx.h>

#include "euclid.hpp"
#include "miller_rabin.hpp"

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

const mpz_class
divmod(
    const mpz_class& Value,
    const mpz_class& Divisor,
    mpz_class& Remainder
);

const mpz_class
modexp(
    const mpz_class& Base,
    const mpz_class& Exponent,
    const mpz_class& Modulus
);

const mpz_class
modexp(
    const mpz_class& Base,
    const int64_t Exponent,
    const mpz_class& Modulus
);

const mpz_class
abs(
    const mpz_class& Value
);

const bool isprime(
    const mpz_class& N
);

const bool isprime(
    const uint64_t N
);

const mpz_class
gcd(
    const mpz_class& A,
    const mpz_class& B
);

const mpz_class
floor_div(
    const mpz_class& A,
    const mpz_class& B
);

const bool
very_fast_prime_test(
    const mpz_class& Candidate
);

const bool
fast_prime_test(
    const mpz_class& Candidate
);

static inline
const bool
divides(
    const mpz_class& N,
    const mpz_class& Value
)
{
    return mpz_divisible_p(N.get_mpz_t(), Value.get_mpz_t());
}

static inline
const bool
divides(
    const mpz_class& N,
    const uint64_t Value
)
{
    return mpz_divisible_ui_p(N.get_mpz_t(), Value);
}

static inline
const bool
divides(
    const __uint128_t N,
    const __uint128_t Value
)
{
    return (N % Value) == 0;
}

static inline
const bool
divides(
    const __uint128_t N,
    const uint64_t Value
)
{
    return (N % Value) == 0;
}

static inline
const bool
divides(
    const uint64_t N,
    const uint64_t Value
)
{
    return (N % Value) == 0;
}

static inline void
increment(
    mpz_class& Value,
    const mpz_class& Amount = 1
)
{
    mpz_add(Value.get_mpz_t(), Value.get_mpz_t(), Amount.get_mpz_t());
}

static inline void
increment(
    mpz_class& Value,
    const uint64_t Amount = 1
)
{
    mpz_add_ui(Value.get_mpz_t(), Value.get_mpz_t(), Amount);
}

static inline void
increment(
    __uint128_t& Value,
    const uint64_t Amount = 1
)
{
    Value += Amount;
}

// Cover all basic types
template <typename T>
static inline void
increment(
    uint64_t& Value,
    const T Amount = 1
)
{
    Value += Amount;
}

static inline void
decrement(
    mpz_class& Value,
    const mpz_class& Amount = 1
)
{
    mpz_sub(Value.get_mpz_t(), Value.get_mpz_t(), Amount.get_mpz_t());
}

static inline void
decrement(
    mpz_class& Value,
    const uint64_t Amount = 1
)
{
    mpz_sub_ui(Value.get_mpz_t(), Value.get_mpz_t(), Amount);
}

static inline void
decrement(
    mpz_class& Value,
    const uint32_t Amount = 1
)
{
    mpz_sub_ui(Value.get_mpz_t(), Value.get_mpz_t(), Amount);
}

// Cover all basic types
template <typename T1, typename T2>
static inline void
decrement(
    T1& Value,
    const T2 Amount = 1
)
{
    Value -= Amount;
}

static inline uint64_t
modulo(
    const mpz_class& A,
    const uint64_t M
)
{
    return mpz_fdiv_ui(A.get_mpz_t(), M);
}

static inline __uint128_t
modulo(
    const __uint128_t& A,
    const uint64_t M
)
{
    return A % M;
}

static inline uint64_t
modulo(
    uint64_t& A,
    const uint64_t M
)
{
    return A % M;
}


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
bit_size(
    const mpz_class& Value
)
{
    return mpz_sizeinbase(Value.get_mpz_t(), 2);
}

template <typename T>
static inline
const size_t
bit_size(
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
    return bit_size(Value) <= 128;
}

static inline __uint128_t
mpz_to_uint128(
    const mpz_class& Value
)
{
    __uint128_t result = Value.get_ui();
    __uint128_t high = mpz_class(Value >> 64).get_ui();
    result |= (high << 64);
    return result;
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

static inline
const bool
IsPerfectSquare(
    const mpz_class& N
)
{
    return mpz_perfect_square_p(N.get_mpz_t());
}

static inline
const bool
IsPerfectSquare(
    const __uint128_t N
)
{
    __uint128_t low = 0;
    __uint128_t high = N;
    while (low <= high) {
        __uint128_t mid = low + (high - low) / 2;
        __uint128_t mid_squared = mid * mid;
        if (mid_squared == N) {
            return true;
        } else if (mid_squared < N) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return false;
}

static inline
const bool
IsPerfectSquare(
    const uint64_t N
)
{
    uint64_t root = static_cast<uint64_t>(std::sqrt(N));
    return root * root == N || (root + 1) * (root + 1) == N;
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