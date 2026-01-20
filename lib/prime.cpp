#include <array>
#include <cstdint>
#include <fstream>
#include <mutex>
#include <span>
#include <vector>

#include "prime.hpp"

namespace primetools {

namespace {
    static const std::array<uint64_t, 100> kSmallPrimes = {
        2, 3, 5, 7, 11, 13, 17, 19, 23, 29,
        31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
        73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
        127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
        179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
        233, 239, 241, 251, 257, 263, 269, 271, 277, 281,
        283, 293, 307, 311, 313, 317, 331, 337, 347, 349,
        353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
        419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
        467, 479, 487, 491, 499, 503, 509, 521, 523, 541
    };

    static std::vector<uint64_t> gCachedPrimesTo;
    static size_t gCachedPrimesToLimit = 0;
    static std::mutex gCachedPrimesMutex;
}

const std::span<const uint64_t>
GetSmallPrimes(
    void
)
{
    return kSmallPrimes;
}

std::span<const uint64_t>
GetPrimesTo(
    const uint64_t Upper
)
{
    std::lock_guard<std::mutex> lock(gCachedPrimesMutex);
    if (Upper <= gCachedPrimesToLimit) {
        // Find the index up to Upper
        size_t index = std::lower_bound(
            gCachedPrimesTo.begin(),
            gCachedPrimesTo.end(),
            Upper + 1
        ) - gCachedPrimesTo.begin();
        return std::span<const uint64_t>(gCachedPrimesTo.data(), index);
    }

    gCachedPrimesTo = GetPrimesInRange<uint64_t>(1, Upper);
    gCachedPrimesToLimit = Upper;
    return std::span<const uint64_t>(gCachedPrimesTo.data(), gCachedPrimesTo.size());
}

const bool LoadPrimeGaps(
    std::string_view FilePath
) {
    std::ifstream ifs(FilePath.data(), std::ios::binary);
    if (!ifs) {
        throw std::runtime_error("Failed to open prime gaps file: " + std::string(FilePath));
    }

    gCachedPrimesTo.clear();

    uint64_t last_prime = 0;
    // Load VLE-encoded gaps
    while (ifs.peek() != EOF) {
        uint64_t gap = 0;
        size_t shift = 0;
        uint8_t byte = 0;
        do {
            ifs.read(reinterpret_cast<char*>(&byte), 1);
            gap |= static_cast<uint64_t>(byte & 0x7F) << shift;
            shift += 7;
        } while (byte & 0x80);
        uint64_t prime = last_prime + gap;
        gCachedPrimesTo.push_back(prime);
        last_prime = prime;
    }

    gCachedPrimesToLimit = gCachedPrimesTo.back();
    return true;
}

} // namespace primetools