#include <array>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <mutex>
#include <span>
#include <thread>
#include <vector>

#include <sys/mman.h>

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
    static uint64_t* gLoadedPrimes = nullptr;
    static size_t gLoadedPrimesCount = 0;
    static std::span<const uint64_t> gLoadedPrimesSpan;
    static size_t gLoadedPrimesLimit = 0;

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

    // Try the loaded primes first
    if (gLoadedPrimes != nullptr && Upper <= gLoadedPrimesLimit) {
        // Do a binary search to find the right size
        size_t left = 0;
        size_t right = gLoadedPrimesCount;
        while (left < right) {
            const size_t mid = left + (right - left) / 2;
            if (gLoadedPrimes[mid] <= Upper) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        return std::span<const uint64_t>(gLoadedPrimes, left);
    }

    if (Upper <= gCachedPrimesToLimit) {
        // Do a binary search to find the right size
        size_t left = 0;
        size_t right = gCachedPrimesTo.size();
        while (left < right) {
            const size_t mid = left + (right - left) / 2;
            if (gCachedPrimesTo[mid] <= Upper) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        return std::span<const uint64_t>(gCachedPrimesTo.data(), left);
    }

    gCachedPrimesTo = GetPrimesInRange<uint64_t>(1, Upper);
    gCachedPrimesToLimit = Upper;
    return gCachedPrimesTo;
}

const bool LoadPrimeGaps(
    std::string_view FilePath
) {
    std::lock_guard<std::mutex> lock(gCachedPrimesMutex);
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

const bool LoadPrimes(
    std::string_view FilePath
) {
    std::lock_guard<std::mutex> lock(gCachedPrimesMutex);

    if (gLoadedPrimes != nullptr) {
        throw std::runtime_error("Primes have already been loaded from a file.");
    }
    
    FILE* fp = fopen(FilePath.data(), "rb");
    if (!fp) {
        throw std::runtime_error("Failed to open primes file: " + std::string(FilePath));
    }

    // Get file size
    size_t file_size = std::filesystem::file_size(FilePath);

    if (gLoadedPrimes == nullptr) {
        size_t num_primes = file_size / sizeof(uint64_t);
        // Memory-map the file
        int fd = fileno(fp);
        void* mapped = mmap(
            nullptr,
            file_size,
            PROT_READ,
            MAP_PRIVATE,
            fd,
            0
        );
        if (mapped == MAP_FAILED) {
            throw std::runtime_error("Failed to memory-map primes file: " + std::string(FilePath));
        }
        gLoadedPrimes = static_cast<uint64_t*>(mapped);
        gLoadedPrimesCount = num_primes;
        gLoadedPrimesSpan = std::span<const uint64_t>(gLoadedPrimes, gLoadedPrimesCount);
        gLoadedPrimesLimit = gLoadedPrimes[num_primes - 1];
    }
    return true;
}

void LoadPrimeGapsInNewThread(
    std::string_view FilePath
) {
    std::thread load_thread([FilePath]() {
        try {
            LoadPrimeGaps(FilePath);
        } catch (const std::exception& e) {
            // Handle error (log it, etc.)
            std::cerr << "Error loading prime gaps: " << e.what() << std::endl;
        }
    });
    load_thread.detach();
}

const size_t
GetCachedPrimesLimit(
    void
) {
    std::lock_guard<std::mutex> lock(gCachedPrimesMutex);
    return gCachedPrimesToLimit;
}

} // namespace primetools