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
#include "prime.hpp"

using namespace primetools;

static const std::string_view HELP_STRING = R"(
Usage: aliquot [options] <number>
Options:
    -p <file>   Load prime gaps from file
    -c <path>   Path to prime factor cache
    -h, --help  Show this help message
)";

static size_t gLastLogLength = 0;
static mpz_class gCurrentTarget;
static size_t gCurrentIndex = 0;

void
HandleLogs(
    std::string_view Message
)
{
    // Erase the last log
    std::cout << '\r' << std::string(gLastLogLength, ' ') << std::flush;
    std::cout << '\r' << std::setw(5) << gCurrentIndex << " : ";
    std::cout << gCurrentTarget << " = " << Message << std::flush;
    gLastLogLength = 5 + 3 + gCurrentTarget.get_str().size() + 3 + Message.size();
}

std::tuple<mpz_class, PrimeFactors<mpz_class>>
SumOfDivisors(
    const mpz_class& N,
    const std::string_view DBPath,
    const size_t NumThreads
)
{
    // Get prime factors of N
    gCurrentTarget = N;
    auto factors = Factorise(N, NumThreads, DBPath, HandleLogs);
    if (!factors) {
        return {0, PrimeFactors<mpz_class>()};
    }

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
    const std::string_view DBPath,
    const bool Verbose,
    const size_t NumThreads
)
{
    // FactorDB<> cache(DBPath);
    std::vector<mpz_class> sequence;
    mpz_class current = N;
    gCurrentIndex = 0;

    while (true) {
        auto [sum, factors] = SumOfDivisors(current, DBPath, NumThreads);
        if (sum == 0) {
            break;
        }
        if (Verbose) {
            HandleLogs(factors.GetString());
            std::cout << std::endl;
        }
        sequence.push_back(sum);
        if (sum == current || DetectLoop(sequence, sum)) {
            break;
        }
        current = sum;
        gCurrentIndex++;
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

    std::string_view db_path, prime_gaps;
    mpz_class number;
    size_t num_threads = 0;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if ((arg == "-p" || arg == "--primes") && i + 1 < argc) {
            prime_gaps = argv[++i];
            // Kick off the loading thread
            primetools::LoadPrimeGapsInNewThread(prime_gaps);
        } else if ((arg == "-d" || arg == "--db") && i + 1 < argc) {
            db_path = argv[++i];
        } else if ((arg == "-t" || arg == "--threads") && i + 1 < argc) {
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
        auto sequence = AliquotSequence(number, db_path, true, num_threads);
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "Error during prime factorization: " << ex.what() << std::endl;
        return 1;
    }

}