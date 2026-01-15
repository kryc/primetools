#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <bit>
#include <cstdint>

namespace primetools
{

template <typename T>
class MiniPRNG
{
public:
    MiniPRNG(T Seed = 0xd34db33f)
    {
        a = 0xf1ea2eed;
        b = c = d = Seed;
        for (int i = 0; i < 20; ++i) {
            Next();
        }
    }
    inline const T Next(void)
    {
        T e;
        if constexpr (std::is_same<T, uint64_t>::value) {
            e = a - std::rotl(b, 7);
            a = b ^ std::rotl(c, 13);
            b = c + std::rotl(d, 37);
        } else if constexpr (std::is_same<T, uint32_t>::value) {
            e = a - std::rotl(b, 32);
            a = b ^ std::rotl(c, 16);
            b = c + std::rotl(d, 11);
        }
        d = e + a;
        c = d + e;
        d = e + a;
        return d;
    }
    inline const T NextEven(void)
    {
        return Next() & ~1ULL;
    }
private:
    T a, b, c, d;
};

// Declare templates specially for 32-bit and 64-bit PRNGs
using MiniPRNG32 = MiniPRNG<uint32_t>;
using MiniPRNG64 = MiniPRNG<uint64_t>;

}

#endif