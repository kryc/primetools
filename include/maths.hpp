#ifndef MATHS_HPP
#define MATHS_HPP

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>

#include <gmpxx.h>

#include "euclid.hpp"

namespace primetools {

static inline __uint128_t
AddMod(
    const __uint128_t A,
    const __uint128_t B,
    const __uint128_t Modulus
)
{
    if (Modulus == 0) {
        throw std::invalid_argument("AddMod: Modulus must be non-zero");
    }

    const __uint128_t a = A % Modulus;
    const __uint128_t b = B % Modulus;

    // Compute (a + b) % Modulus without overflowing __uint128_t
    if (a >= Modulus - b) {
        return a - (Modulus - b);
    }
    return a + b;
}

static inline __uint128_t
MulMod(
    __uint128_t A,
    __uint128_t B,
    const __uint128_t Modulus
)
{
    if (Modulus == 0) {
        throw std::invalid_argument("MulMod: Modulus must be non-zero");
    }

    // Russian peasant multiplication with modular reduction.
    // Avoids needing a 256-bit intermediate for __uint128_t.
    __uint128_t result = 0;
    A %= Modulus;
    while (B > 0) {
        if (B & 1) {
            result = AddMod(result, A, Modulus);
        }
        B >>= 1;
        if (B) {
            A = AddMod(A, A, Modulus);
        }
    }
    return result;
}

const mpz_class
abs(
    const mpz_class& Value
);

// Greatest common divisor (GCD) using Euclidean algorithm
template <typename T>
const T
Gcd(
    const T& A,
    const T& B
)
{
    if constexpr (std::is_same_v<T, mpz_class>) {
        mpz_class result;
        mpz_gcd(result.get_mpz_t(), A.get_mpz_t(), B.get_mpz_t());
        return result;
    } else {
        return EuclideanAlgorithm<T>(A, B);
    }
}

// Integer square root
static inline
const mpz_class
Sqrt(
    const mpz_class& Value
)
{
    mpz_class result;
    mpz_sqrt(result.get_mpz_t(), Value.get_mpz_t());
    return result;
}

static inline
const uint64_t
Sqrt(
    const uint64_t& Value
)
{
    uint64_t x = Value;
    uint64_t y = (x + 1) / 2;
    while (y < x) {
        x = y;
        y = (Value / x + x) / 2;
    }
    return x;
}

static inline
const __uint128_t
Sqrt(
    const __uint128_t& Value
)
{
    __uint128_t x = Value;
    __uint128_t y = (x + 1) / 2;
    while (y < x) {
        x = y;
        y = (Value / x + x) / 2;
    }
    return x;
}

static inline
const mpz_class
Cbrt(
    const mpz_class& Value
)
{
    mpz_class result;
    mpz_root(result.get_mpz_t(), Value.get_mpz_t(), 3);
    return result;
}

static inline
const uint64_t
Cbrt(
    const uint64_t& Value
)
{
    uint64_t x = Value;
    uint64_t y = (2 * x + Value / (x * x)) / 3;
    while (y < x) {
        x = y;
        y = (2 * x + Value / (x * x)) / 3;
    }
    return x;
}

static inline
const __uint128_t
Cbrt(
    const __uint128_t& Value
)
{
    __uint128_t x = Value;
    __uint128_t y = (2 * x + Value / (x * x)) / 3;
    while (y < x) {
        x = y;
        y = (2 * x + Value / (x * x)) / 3;
    }
    return x;
}

static inline
const mpz_class
Root(
    const mpz_class& Value,
    const size_t Root
)
{
    mpz_class result;
    mpz_root(result.get_mpz_t(), Value.get_mpz_t(), Root);
    return result;
}

static inline
const uint64_t
Root(
    const uint64_t& Value,
    const size_t Root
)
{
    if (Root == 0) {
        throw std::invalid_argument("Root: Root must be non-zero");
    }
    if (Root == 1) {
        return Value;
    }
    uint64_t low = 0;
    uint64_t high = Value;
    while (low <= high) {
        uint64_t mid = low + (high - low) / 2;
        uint64_t mid_pow = 1;
        for (size_t i = 0; i < Root; ++i) {
            mid_pow *= mid;
            if (mid_pow > Value) {
                break;
            }
        }
        if (mid_pow == Value) {
            return mid;
        } else if (mid_pow < Value) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return high;
}

static inline
const __uint128_t
Root(
    const __uint128_t& Value,
    const size_t Root
)
{
    if (Root == 0) {
        throw std::invalid_argument("Root: Root must be non-zero");
    }
    if (Root == 1) {
        return Value;
    }
    __uint128_t low = 0;
    __uint128_t high = Value;
    while (low <= high) {
        __uint128_t mid = low + (high - low) / 2;
        __uint128_t mid_pow = 1;
        for (size_t i = 0; i < Root; ++i) {
            mid_pow *= mid;
            if (mid_pow > Value) {
                break;
            }
        }
        if (mid_pow == Value) {
            return mid;
        } else if (mid_pow < Value) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return high;
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
    const uint64_t root = Sqrt(N);
    return root * root == N || (root + 1) * (root + 1) == N;
}

static inline
const bool
IsEven(
    const mpz_class& N
)
{
    return mpz_even_p(N.get_mpz_t());
}

static inline
const bool
IsEven(
    const __uint128_t N
)
{
    return (N & 1) == 0;
}

static inline
const bool
IsEven(
    const uint64_t N
)
{
    return (N & 1) == 0;
}

// Integer exponentiation by squaring
template <typename T>
const T
Pow(
    const T& Base,
    const size_t Exponent
)
{
    if (Exponent == 0) {
        return T(1);
    }
    if constexpr (std::is_same_v<T, mpz_class>) {
        mpz_class result;
        mpz_pow_ui(result.get_mpz_t(), Base.get_mpz_t(), Exponent);
        return result;
    } else {
        T result = 1;
        T base = Base;
        size_t exp = Exponent;

        while (exp > 0) {
            if (exp & 1) {
                if (result > (std::numeric_limits<T>::max() / base)) {
                    throw std::overflow_error("Overflow in Pow function");
                }
                result *= base;
            }
            exp /= 2;
            if (exp > 0) {
                if (base > (std::numeric_limits<T>::max() / base)) {
                    throw std::overflow_error("Overflow in Pow function");
                }
                base *= base;
            }
        }
        return result;
    }
}

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

const mpz_class
floor_div(
    const mpz_class& A,
    const mpz_class& B
);


const mpz_class
divmod(
    const mpz_class& Value,
    const mpz_class& Divisor,
    mpz_class& Remainder
);

const mpz_class
ModExp(
    const mpz_class& Base,
    const mpz_class& Exponent,
    const mpz_class& Modulus
);

const mpz_class
ModExp(
    const mpz_class& Base,
    const int64_t Exponent,
    const mpz_class& Modulus
);

template <typename TE, typename TM>
const __uint128_t
ModExp(
    const __uint128_t& Base,
    const TE& Exponent,
    const TM& Modulus
)
{
    __uint128_t result = 1;
    __uint128_t base = Base % Modulus;
    TE exp = Exponent;

    while (exp > 0) {
        if (exp % 2 == 1) {
            result = MulMod(result, base, static_cast<__uint128_t>(Modulus));
        }
        exp /= 2;
        base = MulMod(base, base, static_cast<__uint128_t>(Modulus));
    }

    return result;
}

} // namespace primetools   

#endif // MATHS_HPP