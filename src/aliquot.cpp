//
// Aliquot
// A simple program to compute aliquot sequences
//

#include <iostream>
#include <iomanip>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <thread>

#include <gmpxx.h>

#include "factorise.hpp"
#include "factors.hpp"

using namespace primetools;

static const std::string_view HELP_STRING = R"(
Usage: aliquot [options] <number>
Options:
    -p <file>   Load prime gaps from file
    -c <path>   Path to prime factor cache
    -h, --help  Show this help message
)";

std::tuple<mpz_class, PrimeFactors<mpz_class>>
SumOfDivisors(
    const mpz_class& N,
    const size_t NumThreads
)
{
    // Get prime factors of N
    auto factors = Factorise(N, NumThreads);
    if (!factors) {
        return {0, PrimeFactors<mpz_class>()};
    }
    // Cache the factors
    // if (Cache.IsOpen()) {
    //     Cache.Write(factors);
    // }
    // Convert the prime factors to a vector of composite factors
    auto composites = factors->GetComposite();
    // Sum the composite factors excluding n itself
    mpz_class sum = 0;
    for (const auto& comp : composites) {
        if (comp != N) {
            sum += comp;
        }
    }
    return {sum, factors.value()};
}

const bool
DetectLoop(
    const std::span<const mpz_class> Sequence,
    const mpz_class& NextValue
)
{
    for (const auto& value : Sequence.subspan(0, Sequence.size() - 1)) {
        if (value == NextValue) {
            return true;
        }
    }
    return false;
}

std::vector<mpz_class>
AliquotSequence(
    const mpz_class& N,
    const std::string_view CachePath,
    const bool Verbose,
    const size_t NumThreads
)
{
    // PrimeFactorCache<> cache(CachePath);
    std::vector<mpz_class> sequence;
    mpz_class current = N;
    size_t index = 0;
    // Output the starting number
    if (Verbose) {
        std::cout <<
            std::setw(5) << index++ <<
            " : " << std::flush;
    }

    while (true) {
        auto [sum, factors] = SumOfDivisors(current, NumThreads);
        if (sum == 0) {
            break;
        }
        if (Verbose) {
            std::cout << factors.GetString() << std::endl;
            std::cout << std::setw(5) << index <<
                " : " <<
                sum << " = " << std::flush;
        }
        sequence.push_back(sum);
        if (sum == current || DetectLoop(sequence, sum)) {
            break;
        }
        current = sum;
        index++;
    }
    if (Verbose) {
        std::cout << std::endl;
    }
    return sequence;
}

int main(
    int argc,
    char* argv[]
) {

    if (argc < 2) {
        std::cerr << HELP_STRING << std::endl;
        return 1;
    }

    std::string_view cache_path;
    mpz_class number;
    size_t num_threads = 0;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        // if ((arg == "-p" || arg == "--primes") && i + 1 < argc) {
        //     prime_gaps = argv[++i];
        //     if (!LoadPrimeGaps(prime_gaps)) {
        //         std::cerr << "Failed to load prime gaps from " << prime_gaps << std::endl;
        //         return 1;
        //     }
        // } else if ((arg == "-c" || arg == "--cache") && i + 1 < argc) {
        //     cache_path = argv[++i];
        if ((arg == "-t" || arg == "--threads") && i + 1 < argc) {
            num_threads = static_cast<size_t>(std::stoul(argv[++i]));
        } else if (arg == "-h" || arg == "--help") {
            std::cout << HELP_STRING << std::endl;
            return 0;
        } else {
            number = argv[i];
        }
    }

    if (number == 0) {
        std::cerr << "Please provide a valid number greater than 0." << std::endl;
        return 1;
    }

    if (num_threads == 0) {
        num_threads = std::thread::hardware_concurrency();
    }

    try {
        std::cout << "Aliquot sequence for " << number << ":" << std::endl;
        auto sequence = AliquotSequence(number, cache_path, true, num_threads);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error during prime factorization: " << ex.what() << std::endl;
        return 1;
    }

}