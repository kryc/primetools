#ifndef TRIAL_DIVISION_HPP
#define TRIAL_DIVISION_HPP

#include <array>
#include <atomic>
#include <future>
#include <iostream>
#include <mutex>
#include <optional>
#include <span>
#include <thread>
#include <variant>
#include <vector>

#include <assert.h>
#include <gmpxx.h>

#include "analyse.hpp"
#include "bigint_avx.hpp"
#include "random.hpp"
#include "trial_division_data.hpp"
#include "util.hpp"

namespace primetools {

const std::span<const __uint128_t>
GetWheelGapsForModulus(
    const size_t Modulus
);

template <typename T>
static inline
const std::tuple<const T, const T, const size_t>
GetUpperAndLowerBounds(
    const T& N,
    const size_t Modulus,
    const bool GuessSize,
    const size_t Bits,
    const T& RangeLower = 0,
    const T& RangeUpper = 0
)
{
    const size_t bits = Bits != 0 ? Bits : (GuessSize ? primetools::GuessSizeOfPrimeFactors(N, true) : primetools::bit_size(N));
    T lower_bound = GuessSize ? (T(1) << (bits - 1)) + RangeLower : RangeLower;
    T upper_bound = (RangeUpper == 0) 
        ? (GuessSize ? primetools::min(T(1) << bits, T(sqrt(N))) : T(sqrt(N))) 
        : (GuessSize ? lower_bound + RangeUpper : RangeUpper);

    // Step back lower_bound to be a multiple of Modulus
    if (lower_bound % Modulus != 0) {
        lower_bound -= (lower_bound % Modulus);
    }

    // Increment upper_bound to be a multiple of Modulus
    if (upper_bound % Modulus != 0) {
        upper_bound += (Modulus - (upper_bound % Modulus));
    }

    return std::make_tuple(lower_bound, upper_bound, bits);
}

template <typename T>
static inline const std::optional<std::pair<T, T>>
AlignCandidateToModulus(
    const T& N,
    T& Candidate,
    const size_t Modulus,
    const bool StepBack = false
)
{
    // Ensure candidate is odd
    if ((Candidate & 1) == 0) {
        Candidate += 1;
    }

    // Align candidate to be congruent to 1 modulo Modulus
    while (primetools::modulo(Candidate, Modulus) != 1 && Candidate > 1) {
        // Check if we found a factor
        if (primetools::divides(N, Candidate) && Candidate > 1) {
            return std::make_pair(Candidate, N / Candidate);
        }
        Candidate += StepBack ? -2 : 2;
    }

    return std::nullopt;
}


template <typename T, const size_t Modulus, const size_t BitSize, typename AT, const size_t Count, const PackingType Packed>
static const std::optional<std::pair<T, T>>
TrialDivisionRange(
    const T& N,
    const T& StartValue,
    const T& EndValue,
    const std::span<const AT, Count> GapArray
)
{
    constexpr size_t WordBits = sizeof(AT) * 8;
    constexpr size_t GapsPerWord =  Packed == FastPack ? (WordBits - 1) / BitSize : WordBits / BitSize;
    constexpr uint64_t GapsMask = Packed == FastPack ? (1 << (BitSize + 1)) - 2 : (1 << BitSize) - 1;

    T candidate = StartValue;

    // Align the candidate to 1 modulo Modulus
    auto align_result = AlignCandidateToModulus<T>(N, candidate, Modulus, true);
    if (align_result.has_value()) {
        return align_result;
    }

    // Special case for wheel as we can use an optimization to rotate
    // the gaps using a bit rotation. This also avoids a branch.
    if constexpr (Modulus == 30)
    {
        uint32_t gapword = WHEEL30GAPSUINT32;
        while (candidate <= EndValue)
        {
            if (primetools::divides(N, candidate) && candidate > 1) {
                return std::make_pair(candidate, N / candidate);
            }

            primetools::increment(candidate, gapword & 0xf);

            gapword = std::rotr(gapword, 4);
        }
    }
    else
    {
        while (candidate <= EndValue)
        {
            for (auto gapword : GapArray)
            {
                for (size_t index = 0; index < GapsPerWord; index++)
                {
                    if (primetools::divides(N, candidate) && candidate > 1) {
                        return std::make_pair(candidate, N / candidate);
                    }

                    if constexpr(Packed == DensePack) {
                        primetools::increment(candidate, (gapword & GapsMask) << 1);
                    } else {
                        primetools::increment(candidate, gapword & GapsMask);
                    }

                    gapword >>= BitSize;
                }
            }
        }
    }

    return std::nullopt;
}

template <typename T, const size_t Modulus, const size_t BitSize, typename AT, const size_t Count, const PackingType Packed>
static const std::optional<std::pair<T, T>>
TrialDivisionRangeSimd(
    const T& N,
    const T& StartValue,
    const T& EndValue,
    const std::span<const AT, Count> GapArray
)
{
#if defined(__AVX512F__)
    constexpr size_t WordBits = sizeof(AT) * 8;
    constexpr size_t GapsPerWord =  Packed == FastPack ? (WordBits - 1) / BitSize : WordBits / BitSize;
    constexpr uint64_t GapsMask = Packed == FastPack ? (1 << (BitSize + 1)) - 2 : (1 << BitSize) - 1;

    std::array<T, 16> candidates;
    candidates[0] = StartValue;

    // Align the candidate to 1 modulo Modulus
    auto align_result = AlignCandidateToModulus<T>(N, candidates[0], Modulus, true);
    if (align_result.has_value()) {
        return align_result;
    }

    std::cout << "Starting SIMD trial division with initial candidate " << candidates[0] << std::endl;

    // Initialize the rest of the candidates
    for (size_t lane = 1; lane < 16; ++lane) {
        candidates[lane] = candidates[lane - 1] + Modulus;
    }

    BigIntAVX<1024> big_candidates;
    big_candidates = candidates;

    if constexpr (Modulus == 30)
    {
        /* TODO*/
    }
    else
    {
        while (true/*candidate <= EndValue*/)
        {
            for (auto gapword : GapArray)
            {
                for (size_t index = 0; index < GapsPerWord; index++)
                {
                    if (big_candidates.divides(N)) {
                        std::cout << "Found divisible candidates!" << std::endl;
                        // Find which lane found the factor
                        for (size_t lane = 0; lane < 16; ++lane) {
                            if (primetools::divides(N, big_candidates.to<T>(lane)) && big_candidates.to<T>(lane) > 1) {
                                T candidate = big_candidates.to<T>(lane);
                                return std::make_pair(candidate, N / candidate);
                            }
                        }
                    }

                    if constexpr(Packed == DensePack) {
                        big_candidates += ((gapword & GapsMask) << 1);
                    } else {
                        big_candidates += (gapword & GapsMask);
                    }

                    gapword >>= BitSize;
                }
            }
            // Increment all candidates by Modulus * Number of lanes (16)
            big_candidates += Modulus * 16;
        }
    }
#endif // __AVX512F__
    return std::nullopt;
}

template <typename T>
static const std::optional<std::pair<T, T>>
TrialDivisionRange(
    const T& N,
    const T& StartValue,
    const T& EndValue,
    const size_t Modulus = 510510
)
{
    switch(Modulus) {
        case 200560490130:
        {
            auto& gaps = GetWheelGapsForModulus(Modulus);
            return TrialDivisionRange<T, 200560490130, 5, __uint128_t, std::dynamic_extent, FastPack>(N, StartValue, EndValue, gaps);
        }
        case 6469693230:
        {
            auto& gaps = GetWheelGapsForModulus(Modulus);
            return TrialDivisionRange<T, 6469693230, 5, __uint128_t, std::dynamic_extent, FastPack>(N, StartValue, EndValue, gaps);
        }
        case 223092870:
            return TrialDivisionRange<T, 223092870, 5, __uint128_t, WHEEL223092870GAP_COUNT, FastPack>(N, StartValue, EndValue, WHEEL223092870GAPS);
        case 9699690:
            return TrialDivisionRange<T, 9699690, 5, __uint128_t, WHEEL9699690GAP_COUNT, FastPack>(N, StartValue, EndValue, WHEEL9699690GAPS);
        case 510510:
            return TrialDivisionRange<T, 510510, 4, uint64_t, WHEEL510510GAP_COUNT, DensePack>(N, StartValue, EndValue, WHEEL510510GAPS);
        case 30030:
            return TrialDivisionRange<T, 30030, 4, uint64_t, WHEEL30030GAP_COUNT, DensePack>(N, StartValue, EndValue, WHEEL30030GAPS);
        case 2310:
            return TrialDivisionRange<T, 2310, 4, uint64_t, WHEEL2310GAP_COUNT, Unpacked>(N, StartValue, EndValue, WHEEL2310GAPS);
        case 210:
            return TrialDivisionRange<T, 210, 4, uint64_t, WHEEL210GAP_COUNT, Unpacked>(N, StartValue, EndValue, WHEEL210GAPS);
        case 30:
            return TrialDivisionRange<T, 30, 4, uint64_t, WHEEL30GAP_COUNT, Unpacked>(N, StartValue, EndValue, WHEEL30GAPS);
        default:
            std::cerr << "Error: Unsupported modulus for trial division: " << Modulus << std::endl;
            return std::nullopt;
    }
}

template <typename T>
static const std::optional<std::pair<T, T>>
TrialDivision(
    const T& N,
    const size_t Modulus,
    const bool GuessSize = true,
    const size_t Bits = 0,
    const T RangeLower = 0,
    const T RangeUpper = 0,
    const bool Simd = false
)
{
    if (N < 2) {
        return std::nullopt;
    }

    // Calculate the bounds and bits
    auto [lower_bound, upper_bound, bits] = GetUpperAndLowerBounds<T>(N, Modulus, GuessSize, Bits, RangeLower, RangeUpper);

    if (GuessSize) {
        std::cout << "Trying factorization of " << bits << "-bit primes using wheel" << Modulus << " factorization." << std::endl;
    }
    std::cout << "Searching primes in range [" << lower_bound << ", " << upper_bound <<
        "] using wheel" << Modulus << " factorization." << std::endl;

    if (Simd) {
        return TrialDivisionRangeSimd<T, 510510, 4, uint64_t, WHEEL510510GAP_COUNT, DensePack>(N, lower_bound, upper_bound, WHEEL510510GAPS);
    } else {
        return TrialDivisionRange<T>(N, lower_bound, upper_bound, Modulus);
    }
}

template <typename T>
const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionBitflip(
    const T& N,
    const size_t MaxIterations,
    const bool UseGuessSize = true
)
{
    if (N < 2) {
        return std::nullopt;
    }

    // Calculate the size of N's constituent primes
    const size_t bits = primetools::GuessSizeOfPrimeFactors(N, true);
    
    // Get the square root of N as this is our upper bound for trial division
    const T upper_bound = T(sqrt(N));
    const size_t upper_bound_bits = primetools::bit_size(upper_bound);

    std::cout << "Trying random factorization of " << bits << "-bit primes " << std::endl;
    std::cout << "WARNING: This is highly inefficient and not recommended!" << std::endl;
    
    // We want to toggle bits between the highest and lowest bit
    const size_t bit_range = UseGuessSize ? bits - 2 : upper_bound_bits - 2; // Exclude highest and lowest bit

    // Make the first candidate the highest and lowest bit set
    T candidate = (T(1) << (bit_range + 1)) | 1;

    // Set up our PRNG
    primetools::MiniPRNG32 prng;

    for (size_t i = 0; i < MaxIterations; ++i) {
        // Flip a random bit in the candidate
        const size_t bit = (prng.Next() % bit_range) + 1;
        primetools::toggle_bit(candidate, bit);

        // Check if it divides N
        if (primetools::divides(N, candidate)) {
            return std::make_pair(candidate, N / candidate);
        }
    }

    return std::nullopt;
};

static constexpr size_t DefaultBlockSize = 1'000'000;

static const size_t
RoundBlockSizeToModulus(
    const size_t BlockSize,
    const size_t Modulus
)
{
    size_t rounded_size = ((BlockSize + Modulus - 1) / Modulus) * Modulus;
    return rounded_size;
}

// Worker function for TrialDivisionMT
template <typename T>
std::optional<std::pair<T, T>>
TrialDivisionMTWorker(
    const T& N,
    const T& LowerBound,
    const T& UpperBound,
    const T& ChunkSize,
    const T& Chunks,
    const size_t ThreadId,
    const size_t NumThreads,
    const size_t Modulus,
    std::atomic<bool>& Found,
    T& CurrentChunk,
    std::mutex& StatusMutex
) {
    T thread_start = LowerBound + (ThreadId * ChunkSize);
    T thread_end = thread_start + ChunkSize;
    while (!Found.load(std::memory_order_relaxed) && thread_start <= UpperBound) {
        T chunk_index = (thread_start - LowerBound) / ChunkSize;
        {
            std::lock_guard<std::mutex> lock(StatusMutex);
            if (chunk_index > CurrentChunk) {
                CurrentChunk = chunk_index;
            }
            std::cout << '\r' << "Chunk " << CurrentChunk << " of " << Chunks << " (" <<
                (CurrentChunk * 100) / Chunks << "%) " << std::flush;
        }
        auto result = TrialDivisionRange<T>(N, thread_start, thread_end, Modulus);
        if (result) {
            Found.store(true, std::memory_order_relaxed);
            return result;
        }
        std::this_thread::yield();
        primetools::increment(thread_start, NumThreads * ChunkSize);
        primetools::increment(thread_end, NumThreads * ChunkSize);
    }
    return std::nullopt;
}

template <typename T>
const std::optional<std::pair<T, T>>
TrialDivisionMT(
    const T& N,
    const size_t Threads,
    const size_t BlockSize,
    const bool GuessSize,
    const size_t Bits,
    const T& RangeLower,
    const T& RangeUpper,
    const size_t Modulus
)
{
    // Use hardware concurrency if Threads == 0
    const size_t num_threads = Threads ? Threads : std::thread::hardware_concurrency();
    const T block_size = RoundBlockSizeToModulus(BlockSize ? BlockSize : DefaultBlockSize, Modulus);

    // Get upper and lower bounds
    auto [lower_bound, upper_bound, bits] = GetUpperAndLowerBounds<T>(N, Modulus, GuessSize, Bits, RangeLower, RangeUpper);

    std::atomic<bool> found{false};
    std::vector<std::future<std::optional<std::pair<T, T>>>> futures;

    // const T chunk_size = Modulus * block_size;
    const T diff = upper_bound - lower_bound;
    const T chunks = (diff / block_size);

    std::cout << "Trying factorization of primes in range [" << primetools::TruncateNumber<T>(lower_bound) << ", " << primetools::TruncateNumber<T>(upper_bound) <<
            "] using modulus " << Modulus << ". " << chunks << " chunks" << std::endl;

    T current_chunk = 0;

    // Mutex for thread-safe status output
    std::mutex status_mutex;

    for (size_t i = 0; i < num_threads; i++) {
        futures.emplace_back(std::async(std::launch::async, TrialDivisionMTWorker<T>,
            std::cref(N),
            std::cref(lower_bound),
            std::cref(upper_bound),
            std::cref(block_size),
            std::cref(chunks),
            i,
            num_threads,
            Modulus,
            std::ref(found),
            std::ref(current_chunk),
            std::ref(status_mutex)
        ));
    }

    std::optional<std::pair<T, T>> result;
    for (auto& fut : futures) {
        auto res = fut.get();
        if (res) {
            result = res;
            break;
        }
    }

    // Terminate the status line
    std::cout << std::endl;
    return result;
}

// Worker function for TrialDivisionRandomMT
template <typename T>
std::optional<std::pair<T, T>>
TrialDivisionRandomMTWorker(
    const T& N,
    const T& LowerBound,
    const T& Chunks,
    const size_t Modulus,
    const T& ChunkSize,
    const size_t ThreadID,
    const size_t NumThreads,
    std::atomic<bool>& Found,
    std::mutex& StatusMutex
) {
    // Split the search space into NumThreads parts
    const T threads_chunks = (Chunks + NumThreads - 1) / NumThreads;
    const T thread_lower = LowerBound + (threads_chunks * ThreadID * ChunkSize);
    const T thread_upper = thread_lower + (threads_chunks * ChunkSize);

    {
        std::lock_guard<std::mutex> lock(StatusMutex);
        std::cout << "Thread " << ThreadID << " searching in range [" << thread_lower << ", " <<
            thread_upper << "] with " << threads_chunks << " chunks." << std::endl;
    }

    // Initialize the PRNG
    MiniPRNG64 prng(ThreadID);

    T start_block = 0;
    while (!Found.load(std::memory_order_relaxed)) {
        primetools::increment(start_block, prng.Next());
        // primetools::increment(start_block, 1);
        if (start_block > threads_chunks) {
            start_block %= threads_chunks;
        }
        T thread_start = thread_lower + (start_block * ChunkSize);
        T thread_end = thread_start + ChunkSize;
        assert(thread_start >= thread_lower && thread_end <= thread_upper);
        // {
        //     std::lock_guard<std::mutex> lock(StatusMutex);
        //     std::cout << '\r' << "Trying block index " << start_block << "/" << threads_chunks << std::flush;
        // }
        auto result = TrialDivisionRange<T>(N, thread_start, thread_end, Modulus);
        if (result) {
            Found.store(true, std::memory_order_relaxed);
            return result;
        }
        std::this_thread::yield();
    }
    return std::nullopt;
}

template <typename T>
const std::optional<std::pair<T, T>>
TrialDivisionRandomMT(
    const T& N,
    const size_t Threads,
    const size_t BlockSize,
    const bool GuessSize,
    const size_t Bits,
    const T& RangeLower,
    const T& RangeUpper,
    const size_t Modulus
)
{
    // Use hardware concurrency if Threads == 0
    const size_t num_threads = Threads ? Threads : std::thread::hardware_concurrency();
    const size_t block_size = RoundBlockSizeToModulus(BlockSize ? BlockSize : DefaultBlockSize, Modulus);

    std::atomic<bool> found{false};
    std::vector<std::future<std::optional<std::pair<T, T>>>> futures;

    // Get bounds and bits
    auto [lower_bound, upper_bound, bits] = GetUpperAndLowerBounds<T>(N, Modulus, GuessSize, Bits, RangeLower, RangeUpper);

    // Calculate the number of Modulus-sized blocks
    T diff = upper_bound - lower_bound;
    T chunks = diff / block_size;

    assert(chunks > 0);

    std::cout << chunks << " chunks of size " << block_size << " (" << block_size << " multiples of Modulus)" << std::endl;

    // Reude thread count if there are fewer chunks than threads
    const size_t effective_threads = chunks < num_threads ? primetools::get_ui(chunks) : num_threads;

    std::mutex status_mutex;

    std::cout << "Trying factorization of " << bits << "-bit primes in range [" << primetools::TruncateNumber(lower_bound) << ", " << primetools::TruncateNumber(upper_bound) <<
            "] using random block search with modulus " << Modulus << ". " << chunks << " chunks" << std::endl;

    for (size_t i = 0; i < effective_threads; i++) {
        futures.emplace_back(std::async(std::launch::async, TrialDivisionRandomMTWorker<T>,
            std::cref(N),
            std::cref(lower_bound),
            std::cref(chunks),
            Modulus,
            block_size,
            i,
            effective_threads,
            std::ref(found),
            std::ref(status_mutex)
        ));
    }

    std::optional<std::pair<T, T>> result;
    for (auto& fut : futures) {
        auto res = fut.get();
        if (res) {
            result = res;
            break;
        }
    }
    
    return result;
}

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionRandom(
    const mpz_class& N,
    const bool GuessSize = true,
    const size_t Bits = 0,
    const mpz_class& RangeLower = 0,
    const mpz_class& RangeUpper = 0,
    const uint64_t Seed = 0,
    const size_t Modulus = 510510,
    const size_t MaxIterations = std::numeric_limits<size_t>::max()
);

const std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionSimd(
    const mpz_class& N,
    const size_t MaxIterations
);

template <typename T>
inline std::optional<std::pair<mpz_class, mpz_class>>
TrialDivisionLinear(const T& N, const size_t Base, const size_t MaxIterations) {
    return TrialDivisionWheel510510<T>(N, MaxIterations);
}

const std::vector<__uint128_t>
GenerateWheelGapsForModulus(
    const size_t Modulus,
    const size_t BitSize = 5,
    const PackingType Packing = PackingType::FastPack
);

} // namespace primetools

#endif // TRIAL_DIVISION_HPP