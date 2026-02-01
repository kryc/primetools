#ifndef TONELLI_SHANKS_HPP
#define TONELLI_SHANKS_HPP

#include <cstdint>
#include <optional>
#include <type_traits>

#include "maths.hpp"

namespace primetools {

// Tonelli–Shanks algorithm for modular square roots.
// Finds x such that x^2 ≡ n (mod p), where p is an odd prime.
// Returns one of the two roots (the other is p-x), or std::nullopt if no root exists.
template <typename T>
static inline std::optional<T>
TonelliShanks(
    T n,
    T p
)
{
    static_assert(std::is_integral_v<T>, "TonelliShanks requires an integral type");
    static_assert(std::is_unsigned_v<T>, "TonelliShanks requires an unsigned type");

    if (p == 2) {
        return static_cast<T>(n & 1u);
    }

    n %= p;
    if (n == 0) {
        return static_cast<T>(0);
    }

    auto modexp = [](T base, T exp, T mod) -> T {
        return static_cast<T>(primetools::ModExp(static_cast<__uint128_t>(base), exp, mod));
    };

    auto modmul = [](T a, T b, T mod) -> T {
        if constexpr (sizeof(T) <= 8) {
            return static_cast<T>((static_cast<__uint128_t>(a) * static_cast<__uint128_t>(b)) % mod);
        } else {
            return static_cast<T>(primetools::MulMod(static_cast<__uint128_t>(a), static_cast<__uint128_t>(b), static_cast<__uint128_t>(mod)));
        }
    };

    // Euler's criterion: n^((p-1)/2) mod p == 1 for quadratic residues.
    const T ls = modexp(n, static_cast<T>((p - 1) / 2), p);
    if (ls != 1) {
        return std::nullopt;
    }

    // Fast path for p ≡ 3 (mod 4).
    if ((p & 3u) == 3u) {
        return modexp(n, static_cast<T>((p + 1) / 4), p);
    }

    // Factor p-1 as q * 2^s with q odd.
    T q = static_cast<T>(p - 1);
    uint32_t s = 0;
    while ((q & 1u) == 0u) {
        q >>= 1u;
        ++s;
    }

    // Find z a quadratic non-residue.
    T z = 2;
    while (z < p) {
        const T zls = modexp(z, static_cast<T>((p - 1) / 2), p);
        if (zls == static_cast<T>(p - 1)) {
            break;
        }
        ++z;
    }
    if (z == p) {
        return std::nullopt;
    }

    T c = modexp(z, q, p);
    T x = modexp(n, static_cast<T>((q + 1) / 2), p);
    T t = modexp(n, q, p);
    uint32_t m = s;

    while (t != 1) {
        uint32_t i = 1;
        T t2i = modmul(t, t, p);
        while (i < m && t2i != 1) {
            t2i = modmul(t2i, t2i, p);
            ++i;
        }
        if (i == m) {
            return std::nullopt;
        }

        T b = c;
        for (uint32_t e = 0; e < (m - i - 1); ++e) {
            b = modmul(b, b, p);
        }

        x = modmul(x, b, p);
        const T b2 = modmul(b, b, p);
        t = modmul(t, b2, p);
        c = b2;
        m = i;
    }

    return x;
}

} // namespace primetools

#endif // TONELLI_SHANKS_HPP
