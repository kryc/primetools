#ifndef TRIAL_DIVISION_HPP
#define TRIAL_DIVISION_HPP

#include <array>
#include <atomic>
#include <functional>
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
#include "factors.hpp"
#include "logging.hpp"
#include "primegenerator.hpp"
#include "random.hpp"
#include "util.hpp"

namespace primetools {

typedef enum TrialDivisionStrategy {
    Linear,
    MeetInTheMiddle,
    Random,
    Hybrid
} TrialDivisionStrategy;

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
    const size_t bits = Bits != 0 ? Bits : (GuessSize ? primetools::GuessSizeOfPrimeFactors(N, true) : primetools::BitSize(N));
    T lower_bound = GuessSize ? (T(1) << (bits - 1)) + RangeLower : RangeLower;
    T upper_bound = (RangeUpper == 0) 
        ? (GuessSize ? primetools::min(T(1) << bits, primetools::Sqrt(N)) : primetools::Sqrt(N)) 
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
static inline void
AlignCandidateToModulus(
    const T& N,
    PrimeFactors<T>& Factors,
    T& Remainder,
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
        if (primetools::divides(N, Candidate) && primetools::IsPrime(Candidate)) {
            while(primetools::divides(Remainder, Candidate)) {
                Factors.AddFactor(Candidate);
                Remainder /= Candidate;
            }
        }
        Candidate += StepBack ? -2 : 2;
    }
}


template <typename T>
static const std::optional<PrimeFactors<T>>
TrialDivisionRange(
    const T& N,
    const T& StartValue,
    const T& EndValue,
    const size_t Modulus
)
{
    T starting_candidate = StartValue;
    T remainder = N;
    PrimeFactors<T> factors;
    // Align the candidate to 1 modulo Modulus
    AlignCandidateToModulus<T>(N, factors, remainder, starting_candidate, Modulus, true);
    // The alignment may have found some factors so check and return
    if (remainder == 1) {
        return factors;
    }

    // Special case for wheel as we can use an optimization to rotate
    // the gaps using a bit rotation. This also avoids a branch.
    if (Modulus == 30)
    {
        uint32_t gapword = kWheel30;
        while (starting_candidate <= EndValue)
        {
            while (primetools::divides(remainder, starting_candidate) && primetools::IsPrime(starting_candidate)) {
                factors.AddFactor(starting_candidate);
                remainder /= starting_candidate;
                if (remainder == 1) {
                    return factors;
                }
            }

            primetools::increment(starting_candidate, gapword & kWheel30Mask);
            gapword = std::rotr(gapword, kWheel30BitsPerGap);
        }
    }
    else
    {
        PossiblePrimeGenerator<T> generator(Modulus, starting_candidate);
        while (generator.Current() <= EndValue)
        {
            const T& candidate = generator.Next();
            if (primetools::divides(N, candidate) && primetools::IsPrime(candidate)) {
                while (primetools::divides(remainder, candidate)) {
                    factors.AddFactor(candidate);
                    remainder /= candidate;
                }
                if (remainder == 1) {
                    break;
                }
            }
        } 
    }

    return factors.Size() > 0 ? std::optional<PrimeFactors<T>>(factors) : std::nullopt;
}

// template <typename T, const size_t Modulus, const size_t BitSize, typename AT, const size_t Count, const PackingType Packed>
// static const std::optional<PrimeFactors<T>>
// TrialDivisionRangeSimd(
//     const T& N,
//     const T& StartValue,
//     const T& EndValue,
//     const std::span<const AT, Count> GapArray
// )
// {
// #if defined(__AVX512F__)
//     constexpr size_t WordBits = sizeof(AT) * 8;
//     constexpr size_t GapsPerWord =  Packed == FastPack ? (WordBits - 1) / BitSize : WordBits / BitSize;
//     constexpr uint64_t GapsMask = Packed == FastPack ? (1 << (BitSize + 1)) - 2 : (1 << BitSize) - 1;

//     std::array<T, 16> candidates;
//     candidates[0] = StartValue;

//     // Align the candidate to 1 modulo Modulus
//     auto align_result = AlignCandidateToModulus<T>(N, candidates[0], Modulus, true);
//     if (align_result.has_value()) {
//         return align_result;
//     }

//     std::cout << "Starting SIMD trial division with initial candidate " << candidates[0] << std::endl;

//     // Initialize the rest of the candidates
//     for (size_t lane = 1; lane < 16; ++lane) {
//         candidates[lane] = candidates[lane - 1] + Modulus;
//     }

//     BigIntAVX<1024> big_candidates;
//     big_candidates = candidates;

//     if constexpr (Modulus == 30)
//     {
//         /* TODO*/
//     }
//     else
//     {
//         while (true/*candidate <= EndValue*/)
//         {
//             for (auto gapword : GapArray)
//             {
//                 for (size_t index = 0; index < GapsPerWord; index++)
//                 {
//                     if (big_candidates.divides(N)) {
//                         std::cout << "Found divisible candidates!" << std::endl;
//                         // Find which lane found the factor
//                         for (size_t lane = 0; lane < 16; ++lane) {
//                             if (primetools::divides(N, big_candidates.to<T>(lane)) && big_candidates.to<T>(lane) > 1) {
//                                 T candidate = big_candidates.to<T>(lane);
//                                 return std::make_pair(candidate, N / candidate);
//                             }
//                         }
//                     }

//                     if constexpr(Packed == DensePack) {
//                         big_candidates += ((gapword & GapsMask) << 1);
//                     } else {
//                         big_candidates += (gapword & GapsMask);
//                     }

//                     gapword >>= BitSize;
//                 }
//             }
//             // Increment all candidates by Modulus * Number of lanes (16)
//             big_candidates += Modulus * 16;
//         }
//     }
// #endif // __AVX512F__
//     return std::nullopt;
// }

template <typename T>
static const std::optional<PrimeFactors<T>>
TrialDivisionLinear(
    const T& N,
    const size_t Modulus,
    const bool GuessSize = true,
    const size_t Bits = 0,
    const T RangeLower = 0,
    const T RangeUpper = 0,
    const bool Simd = false,
    const bool Status = false
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

    if (Status) {
        std::cout << "Searching primes in range [" << lower_bound << ", " << upper_bound <<
            "] using wheel" << Modulus << " factorization." << std::endl;
    }

    // if (Simd) {
    //     return TrialDivisionRangeSimd<T, 510510, 4, uint64_t, WHEEL510510GAP_COUNT, DensePack>(N, lower_bound, upper_bound, WHEEL510510GAPS);
    // } else {
        return TrialDivisionRange<T>(N, lower_bound, upper_bound, Modulus);
    // }
}

// template <typename T>
// const std::optional<PrimeFactors<T>>
// TrialDivisionBitflip(
//     const T& N,
//     const size_t MaxIterations,
//     const bool UseGuessSize = true
// )
// {
//     if (N < 2) {
//         return std::nullopt;
//     }

//     // Calculate the size of N's constituent primes
//     const size_t bits = primetools::GuessSizeOfPrimeFactors(N, true);
    
//     // Get the square root of N as this is our upper bound for trial division
//     const T upper_bound = T(sqrt(N));
//     const size_t upper_bound_bits = primetools::BitSize(upper_bound);

//     std::cout << "Trying random factorization of " << bits << "-bit primes " << std::endl;
//     std::cout << "WARNING: This is highly inefficient and not recommended!" << std::endl;
    
//     // We want to toggle bits between the highest and lowest bit
//     const size_t bit_range = UseGuessSize ? bits - 2 : upper_bound_bits - 2; // Exclude highest and lowest bit

//     // Make the first candidate the highest and lowest bit set
//     T candidate = (T(1) << (bit_range + 1)) | 1;

//     // Set up our PRNG
//     primetools::MiniPRNG32 prng;

//     for (size_t i = 0; i < MaxIterations; ++i) {
//         // Flip a random bit in the candidate
//         const size_t bit = (prng.Next() % bit_range) + 1;
//         primetools::toggle_bit(candidate, bit);

//         // Check if it divides N
//         if (primetools::divides(N, candidate)) {
//             return PrimeFactors<T>::FromPair(candidate, N / candidate);
//         }
//     }

//     return std::nullopt;
// };

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

template <typename T>
const size_t
TrialDivisionLinearWorker(
    const T& N,
    PrimeFactors<T>& Factors,
    const T& LowerBound,
    const T& UpperBound,
    const T& ChunkSize,
    const T& Chunks,
    const size_t ThreadId,
    const size_t NumThreads,
    const size_t Modulus,
    std::atomic<bool>& Found,
    T& CurrentChunk,
    std::mutex& StatusMutex,
    const LogCallback LogFn = primetools::LogQuiet
) {
    T thread_start = LowerBound + (ThreadId * ChunkSize);
    T thread_end = thread_start + ChunkSize;
    size_t count = 0;

    // Make a thread-local copy of Factors and Remainder to avoid contention
    PrimeFactors<T> thread_factors;
    T thread_remainder;

    while (!Found.load() && thread_start <= UpperBound) {
        {
            std::lock_guard<std::mutex> lock(StatusMutex);
            thread_remainder = N / Factors.Product();
        }
        // Check if remainder is prime. This can happen as numbers have
        // At most one prime factor greater than sqrt(N) and we usually
        // only do trial division up to sqrt(N). We can return early if
        // we have found all other factors and only the prime remainder is left.
        // Also check its above the upper bound to avoid adding twice.
        if (primetools::IsPrime(thread_remainder) && thread_remainder > UpperBound) {
            std::lock_guard<std::mutex> lock(StatusMutex);
            if (!Factors.HasFactor(thread_remainder)) {
                Factors.AddFactor(thread_remainder);
                count += 1;
            }
            Found.store(true);
            break;
        }
        if (thread_remainder == 1) {
            Found.store(true);
            break;
        }

        // Report our status
        T chunk_index = (thread_start - LowerBound) / ChunkSize;
        {
            std::lock_guard<std::mutex> lock(StatusMutex);
            if (chunk_index > CurrentChunk) {
                CurrentChunk = chunk_index;
            }
            LogFn("Chunk " + primetools::ToString(CurrentChunk) +
            " (" + primetools::ToString(thread_start) + " to " + primetools::ToString(thread_end) + ")" +
            " of " + primetools::ToString(Chunks) + " (" +
                primetools::ToString((CurrentChunk * 100) / Chunks) + "%)");
        }
        
        auto result = TrialDivisionRange<T>(thread_remainder, thread_start, thread_end, Modulus);
        if (result) {
            // We found some factors. Merge them back into the global Factors
            // and update our local remainder
            std::lock_guard<std::mutex> lock(StatusMutex);
            count += result->Count();
            Factors.Update(result.value());
            if (Factors.Product() == N) {
                Found.store(true);
            }
        }
        std::this_thread::yield();
        primetools::increment(thread_start, NumThreads * ChunkSize);
        primetools::increment(thread_end, NumThreads * ChunkSize);
    }
    return count;
}

// Worker function for TrialDivisionRandomMT
template <typename T>
const size_t
TrialDivisionRandomWorker(
    const T& N,
    PrimeFactors<T>& Factors,
    const T& LowerBound,
    const T& UpperBound,
    const T& ChunkSize,
    const T& Chunks,
    const size_t ThreadId,
    const size_t NumThreads,
    const size_t Modulus,
    std::atomic<bool>& Found,
    T& CurrentChunk,
    std::mutex& StatusMutex,
    const LogCallback LogFn = primetools::LogQuiet
) {
    // Make a thread-local copy of Factors and Remainder to avoid contention
    PrimeFactors<T> thread_factors;
    T thread_remainder;
    size_t count = 0;

    {
        std::lock_guard<std::mutex> lock(StatusMutex);
        thread_remainder = N / Factors.Product();
    }
    if (thread_remainder == 1) {
        Found.store(true);
        return 0;
    }

    // Initialize the PRNG
    MiniPRNG64 prng(ThreadId);

    T chunk = 0;

    while (!Found.load()) {
        {
            std::lock_guard<std::mutex> lock(StatusMutex);
            thread_remainder = N / Factors.Product();
        }
        // Check if remainder is prime. This can happen as numbers have
        // At most one prime factor greater than sqrt(N) and we usually
        // only do trial division up to sqrt(N). We can return early if
        // we have found all other factors and only the prime remainder is left.
        if (primetools::IsPrime(thread_remainder) && thread_remainder > UpperBound) {
            std::lock_guard<std::mutex> lock(StatusMutex);
            if (!Factors.HasFactor(thread_remainder)) {
                Factors.AddFactor(thread_remainder);
                count += 1;
            }
            Found.store(true);
            break;
        }
        if (thread_remainder == 1) {
            Found.store(true);
            break;
        }

        const uint64_t chunkOffset = prng.Next();
        primetools::increment(chunk, chunkOffset * NumThreads);
        if (chunk > Chunks) {
            chunk %= Chunks;
        }
        T thread_start = LowerBound + (chunk * ChunkSize);
        T thread_end = thread_start + ChunkSize;
        assert(thread_start >= LowerBound && thread_end <= UpperBound);
        
        auto result = TrialDivisionRange<T>(thread_remainder, thread_start, thread_end, Modulus);
        if (result) {
            // We found some factors. Merge them back into the global Factors
            // and update our local remainder
            std::lock_guard<std::mutex> lock(StatusMutex);
            count += result->Count();
            Factors.Update(result.value());
            if (Factors.Product() == N) {
                Found.store(true);
            }
        }
        std::this_thread::yield();
    }

    return count;
}

template <typename T>
const size_t
TrialDivisionMeetInTheMiddleWorker(
    const T& N,
    PrimeFactors<T>& Factors,
    const T& LowerBound,
    const T& UpperBound,
    const T& ChunkSize,
    const T& Chunks,
    const size_t ThreadId,
    const size_t NumThreads,
    const size_t Modulus,
    std::atomic<bool>& Found,
    T& CurrentChunk,
    std::mutex& StatusMutex,
    const LogCallback LogFn = primetools::LogQuiet
) {
    // Meet in the middle works by threads starting at either the lower bound
    // and working up, or starting at the upper bound and working down.
    // Even thread IDs start start from the lower bound, odd thread IDs start
    // from the upper bound.

    const size_t lower_threads = (NumThreads + 1) / 2; // even ThreadId values
    const size_t upper_threads = NumThreads / 2;       // odd ThreadId values

    T thread_start;
    if (ThreadId % 2 == 0) {
        thread_start = LowerBound + T(ThreadId / 2) * ChunkSize;
    } else {
        thread_start = UpperBound - T((ThreadId / 2) + 1) * ChunkSize;
    }
    T thread_end = primetools::min(thread_start + ChunkSize, UpperBound);

    // Work out the middle point to avoid overlapping ranges
    // Once thread_start exceeds middle_point, we can stop
    const T middle_point = LowerBound + ((UpperBound - LowerBound) / 2);
    
    size_t count = 0;

    // Make a thread-local copy of Factors and Remainder to avoid contention
    PrimeFactors<T> thread_factors;
    T thread_remainder;

    while (
        !Found.load() && 
        ((ThreadId % 2 == 0 && thread_start <= middle_point) ||
         (ThreadId % 2 == 1 && thread_start >= middle_point)))
    {
        {
            std::lock_guard<std::mutex> lock(StatusMutex);
            thread_remainder = N / Factors.Product();
        }
        // Check if remainder is prime. This can happen as numbers have
        // At most one prime factor greater than sqrt(N) and we usually
        // only do trial division up to sqrt(N). We can return early if
        // we have found all other factors and only the prime remainder is left.
        if (primetools::IsPrime(thread_remainder) && thread_remainder > UpperBound) {
            std::lock_guard<std::mutex> lock(StatusMutex);
            if (!Factors.HasFactor(thread_remainder)) {
                Factors.AddFactor(thread_remainder);
                count += 1;
            }
            Found.store(true);
            break;
        }
        if (thread_remainder == 1) {
            Found.store(true);
            break;
        }
        
        // Report our status if thread 0
        if (ThreadId == 0) {
            std::lock_guard<std::mutex> lock(StatusMutex);
            thread_remainder = N / Factors.Product();
            T chunk_index = (thread_start - LowerBound) / ChunkSize;
            
            if (chunk_index > CurrentChunk) {
                CurrentChunk = chunk_index;
            }
            
            // Thread zero only reports up to the middle
            T thread0_chunks = Chunks / 2;
            T percent_complete = (CurrentChunk * 100) / thread0_chunks;
            LogFn(primetools::ToString(CurrentChunk) + " / " + primetools::ToString(Chunks / 2) + " (" + primetools::ToString(percent_complete) + "%)");
        }

        auto result = TrialDivisionRange<T>(thread_remainder, thread_start, thread_end, Modulus);
        if (result) {
            // We found some factors. Merge them back into the global Factors
            // and update our local remainder
            std::lock_guard<std::mutex> lock(StatusMutex);
            count += result->Count();
            Factors.Update(result.value());
            if (Factors.Product() == N) {
                Found.store(true);
            }
        }
        std::this_thread::yield();
        if (ThreadId % 2 == 0) {
            primetools::increment(thread_start, lower_threads * ChunkSize);
            thread_end = primetools::min(thread_start + ChunkSize, UpperBound);
        } else {
            // upper_threads is non-zero for any odd ThreadId that exists.
            primetools::decrement(thread_start, upper_threads * ChunkSize);
            thread_end = primetools::min(thread_start + ChunkSize, UpperBound);
        }
    }
    return count;
}

template <typename T>
std::optional<PrimeFactors<T>>
TrialDivision(
    const T& N,
    const size_t Threads,
    const size_t BlockSize,
    const bool GuessSize,
    const size_t Bits,
    const T& RangeLower,
    const T& RangeUpper,
    const size_t Modulus,
    const TrialDivisionStrategy Strategy = TrialDivisionStrategy::Linear,
    const LogCallback LogFn = primetools::LogQuiet
)
{
    const size_t num_threads = Threads ? Threads : std::thread::hardware_concurrency();
    const T block_size = RoundBlockSizeToModulus(BlockSize ? BlockSize : DefaultBlockSize, Modulus);

    // Get upper and lower bounds
    auto [lower_bound, upper_bound, bits] = GetUpperAndLowerBounds<T>(N, Modulus, GuessSize, Bits, RangeLower, RangeUpper);

    if (upper_bound < lower_bound) {
        throw std::runtime_error("Trial division upper bound is less than lower bound.");
    }

    PrimeFactors<T> factors;

    std::atomic<bool> found{false};
    std::vector<std::future<size_t>> futures;

    const T diff = upper_bound - lower_bound;
    // Number of chunk start positions where start <= upper_bound.
    const T chunks = (diff / block_size) + 1;

    LogFn("Trying factorization of primes in range [" + primetools::TruncateNumber<T>(lower_bound) + ", " + primetools::TruncateNumber<T>(upper_bound) +
        "] using modulus " + std::to_string(Modulus) + ". " + primetools::ToString(chunks) + " chunks");

    // Single-threaded linear search case
    if (num_threads == 1 && Strategy == TrialDivisionStrategy::Linear) {
        auto result = TrialDivisionRange<T>(N, lower_bound, upper_bound, Modulus);
        if (result) {
            return result;
        } else {
            return std::nullopt;
        }
    }

    T current_chunk = 0;

    // Mutex for thread-safe status output
    std::mutex status_mutex;

    if (Strategy == TrialDivisionStrategy::Linear) {
        for (size_t i = 0; i < num_threads; i++) {
            futures.emplace_back(std::async(std::launch::async, TrialDivisionLinearWorker<T>,
                std::ref(N), std::ref(factors),
                std::cref(lower_bound), std::cref(upper_bound),
                std::cref(block_size), std::cref(chunks),
                i,
                num_threads,
                Modulus,
                std::ref(found),
                std::ref(current_chunk),
                std::ref(status_mutex),
                LogFn
            ));
        }
    } else if (Strategy == TrialDivisionStrategy::MeetInTheMiddle) {
        for (size_t i = 0; i < num_threads; i++) {
            futures.emplace_back(std::async(std::launch::async, TrialDivisionMeetInTheMiddleWorker<T>,
                std::ref(N), std::ref(factors),
                std::cref(lower_bound), std::cref(upper_bound),
                std::cref(block_size), std::cref(chunks),
                i,
                num_threads,
                Modulus,
                std::ref(found),
                std::ref(current_chunk),
                std::ref(status_mutex),
                LogFn
            ));
        }
    } else if (Strategy == TrialDivisionStrategy::Random) {
        for (size_t i = 0; i < num_threads; i++) {
            futures.emplace_back(std::async(std::launch::async, TrialDivisionRandomWorker<T>,
                std::ref(N), std::ref(factors),
                std::cref(lower_bound), std::cref(upper_bound),
                std::cref(block_size), std::cref(chunks),
                i,
                num_threads,
                Modulus,
                std::ref(found),
                std::ref(current_chunk),
                std::ref(status_mutex),
                LogFn
            ));
        }
    } else if (Strategy == TrialDivisionStrategy::Hybrid) {
        // We divide the number of threads in half for linear and half for random
        // Use the larger of the two for linear to ensure coverage
        const size_t linear_threads = (num_threads + 1) / 2;
        const size_t random_threads = num_threads - linear_threads;
        // If we have at least two linear threads, we can use meet-in-the-middle
        if (linear_threads >= 2) {
            for (size_t i = 0; i < linear_threads; i++) {
                futures.emplace_back(std::async(std::launch::async, TrialDivisionMeetInTheMiddleWorker<T>,
                    std::ref(N), std::ref(factors),
                    std::cref(lower_bound), std::cref(upper_bound),
                    std::cref(block_size), std::cref(chunks),
                    i,
                    linear_threads,
                    Modulus,
                    std::ref(found),
                    std::ref(current_chunk),
                    std::ref(status_mutex),
                    LogFn
                ));
            }
        } else {
            for (size_t i = 0; i < linear_threads; i++) {
                futures.emplace_back(std::async(std::launch::async, TrialDivisionLinearWorker<T>,
                    std::ref(N), std::ref(factors),
                    std::cref(lower_bound), std::cref(upper_bound),
                    std::cref(block_size), std::cref(chunks),
                    i,
                    linear_threads,
                    Modulus,
                    std::ref(found),
                    std::ref(current_chunk),
                    std::ref(status_mutex),
                    LogFn
                ));
            }
        }
        for (size_t i = 0; i < random_threads; i++) {
            futures.emplace_back(std::async(std::launch::async, TrialDivisionRandomWorker<T>,
                std::ref(N), std::ref(factors),
                std::cref(lower_bound), std::cref(upper_bound),
                std::cref(block_size), std::cref(chunks),
                i,
                random_threads,
                Modulus,
                std::ref(found),
                std::ref(current_chunk),
                std::ref(status_mutex),
                LogFn
            ));
        }
    }

    // size_t result = 0;
    for (auto& fut : futures) {
        fut.get();
    }

    // Terminate the status line
    return factors.Size() > 0 ? std::optional<PrimeFactors<T>>(factors) : std::nullopt;
}

} // namespace primetools

#endif // TRIAL_DIVISION_HPP