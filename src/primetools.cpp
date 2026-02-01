#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <string_view>

#include <stdint.h>
#include <gmpxx.h>

#include "factorise.hpp"
#include "factors.hpp"
#include "prime.hpp"
#include "util.hpp"

template <typename T>
void
OutputFactors(
    const std::optional<primetools::PrimeFactors<T>>& Factors
)
{
    if (Factors) {
        std::cout << Factors->GetString() << std::endl;
    } else {
        std::cout << "No factors found" << std::endl;
    }
}

template <typename T>
void
OutputFactors(
    const std::optional<std::pair<T, T>>& Factors
)
{
    if (Factors) {
        std::cout << Factors->first << " * " << Factors->second << std::endl;
    } else {
        std::cout << "No factors found" << std::endl;
    }
}

int main(
    int argc,
    char* argv[]
)
{
    // Check if the number of arguments is correct
    if (argc < 2) {
        return 1;
    }

    std::map<std::string_view, std::string_view> flags;
    std::vector<std::string_view> positionals;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (arg.rfind("--", 0) == 0) {
            size_t eq_pos = arg.find('=');
            if (eq_pos != std::string::npos) {
                std::string_view key = arg.substr(2, eq_pos - 2);
                std::string_view value = arg.substr(eq_pos + 1);
                flags[key] = value;
            } else {
                std::string_view key = arg.substr(2);
                flags[key] = "true";
            }
        } else {
            positionals.push_back(arg);
        }
    }

    if (positionals.size() < 1) {
        std::cerr << "No action specified." << std::endl;
        return 1;
    }

    const std::string_view action = positionals[0];

    if (action == "factorise" || action == "factorize")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " factorise <number>" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const bool quiet = flags.find("quiet") != flags.end();
        const std::string_view factor_db_path = flags.find("db") != flags.end() ? flags["db"] : "";

        const mpz_class n(positionals[1].data());

        auto result = primetools::Factorise(n, 0, factor_db_path, quiet ? primetools::LogQuiet : primetools::LogStdOut);
        OutputFactors(result);
    }
    else if (action == "quadratic")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " quadratic <number>" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const bool quiet = flags.find("quiet") != flags.end();
        const mpz_class n(positionals[1].data());
        std::cout << "Quadratic Sieve factorization of " << primetools::TruncateNumber(n) << std::endl;
        auto result = primetools::QuadraticSieveFactor(n, quiet ? primetools::LogQuiet : primetools::LogStdOut);
        OutputFactors(result);
    }
    else if (action == "fermat")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " fermat <number> [--fermat|--ffa|--mffv2|--mffv4] [--max-iterations=<max>]" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const mpz_class n(positionals[1].data());
        const size_t max = flags.find("max-iterations") != flags.end() ? std::stoul(std::string(flags["max-iterations"])) : std::numeric_limits<size_t>::max();
        const std::string_view algorithm_str = flags.find("algorithm") != flags.end() ? flags["algorithm"] : "fermat";

        // Parse the algorithm flag
        const primetools::FermatAlgorithm algorithm = primetools::GetFermatAlgorithmFromString(algorithm_str);

        std::cout << "Using " << primetools::FermatAlgorithmToString(algorithm) << " factorization of " << primetools::TruncateNumber(n) << std::endl;
        auto result = primetools::FermatFactorisation(n, 0, algorithm, 0, max);
        OutputFactors(result);
    }
    else if (action == "rho")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " rho <number> [--pollard]" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const mpz_class n(positionals[1].data());
        const bool pollard = flags.find("pollard") != flags.end();
        const size_t threads = flags.find("threads") != flags.end() ? std::stoul(std::string(flags["threads"])) : 1;

        if (pollard) {
            std::cout << "Pollard's rho factorization of " << primetools::TruncateNumber(n) << std::endl;
            auto result = primetools::PollardsRho(n);
            OutputFactors(result);
        }
        else {
            std::cout << "Brent-Pollard's rho factorization of " << primetools::TruncateNumber(n) << std::endl;
            if (threads > 1) {
                auto result = primetools::BrentPollardsRhoMT(n, threads);
                OutputFactors(result);    
            } else {
                auto result = primetools::BrentPollardsRho(n);
                OutputFactors(result);
            }
        }
    }
    else if (action == "pminus1" || action == "p-1" || action == "pollardp1" || action == "p1" || action == "pollardp-1")
    {
        if (argc < 2) {
            std::cerr << "Usage: " << argv[0] << " pminus1 <number> [--bound=<B>] [--bases=<bases>]" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const mpz_class n(positionals[1].data());
        const size_t bound = flags.find("bound") != flags.end() ? std::stoul(std::string(flags["bound"])) : std::powl(2, 32);
        const size_t bases = flags.find("bases") != flags.end() ? std::stoul(std::string(flags["bases"])) : 1'000'000;
        auto result = primetools::PollardsPMinus1(n, bound, bases);
        OutputFactors(result);
    }
    else if (action == "squfof")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " squfof <number>" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const mpz_class n(positionals[1].data());
        auto result = primetools::SQUFOF(n);
        OutputFactors(result);
    }
    else if (action == "random")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " random <number> [max_iterations]" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        // Check if the user specified a modulus
        size_t modulus = 510510; // Default to 510510
        if (flags.find("modulus") != flags.end() || flags.find("wheel") != flags.end()) {
            modulus = flags.find("modulus") != flags.end() ? std::stoul(std::string(flags["modulus"])) : std::stoul(std::string(flags["wheel"]));
        }

        const mpz_class n(positionals[1].data());
        const mpz_class start = flags.find("start") != flags.end() ? mpz_class(flags["start"].data()) : mpz_class(1);
        const mpz_class end = flags.find("end") != flags.end() ? mpz_class(flags["end"].data()) : mpz_class(0);
        const bool noguess = flags.find("no-guess") != flags.end();
        const size_t bits = flags.find("bits") != flags.end() ? std::stoul(std::string(flags["bits"])) : 0;
        const size_t threads = flags.find("threads") != flags.end() ? std::stoul(std::string(flags["threads"])) : 1;
        const size_t blocksize = flags.find("blocksize") != flags.end() ? std::stoul(std::string(flags["blocksize"])) : 0;
        // const uint64_t seed = flags.find("seed") != flags.end() ? std::stoull(std::string(flags["seed"])) : 0;
        // const size_t max_iterations = flags.find("max-iterations") != flags.end() ? std::stoul(std::string(flags["max-iterations"])) : std::numeric_limits<size_t>::max();

        std::cout << "Random prime factorization of " << primetools::TruncateNumber(n) << std::endl;

        auto result = primetools::TrialDivision<mpz_class>(n, threads, blocksize, !noguess, bits, start, end, modulus, primetools::TrialDivisionStrategy::Random);
        OutputFactors(result);
    }
    else if (action == "trial" || action == "wheel" || action == "trialdivision" || action == "td")
    {
        // Check and get the number to factor
        if (positionals.size() < 2) {
            std::cerr << "Usage: " << argv[0] << " trial <number> [--modulus=<modulus>] [--no-guess] [--bits=<bits>] [--start=<start>] [--end=<end>] [--simd] [--use-simd] [--max-iterations=<max>]" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        // Check if the user specified a modulus
        size_t modulus = 510510; // Default to 510510
        if (flags.find("modulus") != flags.end() || flags.find("wheel") != flags.end()) {
            modulus = flags.find("modulus") != flags.end() ? std::stoul(std::string(flags["modulus"])) : std::stoul(std::string(flags["wheel"]));
        }

        const mpz_class n(positionals[1].data());
        const mpz_class start = flags.find("start") != flags.end() ? mpz_class(flags["start"].data()) : mpz_class(0);
        const mpz_class end = flags.find("end") != flags.end() ? mpz_class(flags["end"].data()) : mpz_class(0);
        const bool noguess = flags.find("no-guess") != flags.end();
        const size_t bits = flags.find("bits") != flags.end() ? std::stoul(std::string(flags["bits"])) : 0;
        // const bool simd = flags.find("simd") != flags.end() || flags.find("use-simd") != flags.end();
        const size_t threads = flags.find("threads") != flags.end() ? std::stoul(std::string(flags["threads"])) : 1;
        const size_t blocksize = flags.find("blocksize") != flags.end() ? std::stoul(std::string(flags["blocksize"])) : 0;

        /*if (simd) {
            const size_t max_iterations = flags.find("max-iterations") != flags.end() ? std::stoul(flags["max-iterations"]) : std::numeric_limits<size_t>::max();
            auto result = primetools::TrialDivisionSimd(n, max_iterations);
            OutputFactors(result);
        }
        else */
        if (primetools::fits_uint128(n)) {
            auto result = primetools::TrialDivision<__uint128_t>(primetools::MpzToUint128(n), threads, blocksize, !noguess, bits, primetools::MpzToUint128(start), primetools::MpzToUint128(end), modulus);
            OutputFactors(result);
        }
        else
        {
            auto result = primetools::TrialDivision(n, threads, blocksize, !noguess, bits, start, end, modulus, primetools::TrialDivisionStrategy::MeetInTheMiddle, primetools::LogStdOut);
            OutputFactors(result);
        }
    }
    // else if (action == "bitflip")
    // {
    //     if (argc < 3) {
    //         std::cerr << "Usage: " << argv[0] << " trialdivision <number> [max_iterations]" << std::endl;
    //         return 1;
    //     }

    //     if (!primetools::is_numeric(positionals[1])) {
    //         std::cerr << "Error: Input is not a valid number." << std::endl;
    //         return 1;
    //     }

    //     const mpz_class n(positionals[1]);
    //     const size_t max_iterations = argc == 5 ? std::stoul(argv[4]) : std::numeric_limits<size_t>::max();

    //     auto result = primetools::TrialDivisionBitflip(n, max_iterations);
    //     OutputFactors(result);
    // }
    else if (action == "isprime")
    {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " isprime <number>" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const mpz_class n(positionals[1].data());

        bool is_prime = primetools::IsPrime(n);
        std::cout << primetools::TruncateNumber(n) << (is_prime ? " is prime." : " is not prime.") << std::endl;
    }
    else if (action == "calculatefermatiterations")
    {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " calculatefermatiterations <number>" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const mpz_class n(positionals[1].data());
        size_t iterations = primetools::CalculateFermatIterations(n);
        std::cout << "Fermat iterations for " << primetools::TruncateNumber(n) << ": " << iterations << std::endl;
    }
    else if (action == "moduli")
    {
        std::cout << "Supported wheel moduli: 30, 210, 2310, 30030, 510510, 9699690, 223092870, 6469693230, 200560490130, 7420738134810, 304250263527210" << std::endl;
    }
    else if (action == "getmodulus")
    {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " getmodulus <prime>" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const uint64_t prime = std::stoull(positionals[1].data());
        if (prime < 2) {
            std::cerr << "Error: Prime must be greater than 1." << std::endl;
            return 1;
        }

        uint64_t moduli = primetools::GetTrialDivisionModuliForPrime(prime);
        std::cout << "Trial division moduli for prime " << prime << ": " << moduli << std::endl;
    }
    else if (action == "nthprime")
    {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " nthprime <n>" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const size_t n = std::stoul(positionals[1].data());
        const mpz_class prime = primetools::GetNthPrime(n);
        std::cout << "The " << n << (n % 10 == 1 && n % 100 != 11 ? "st" : (n % 10 == 2 && n % 100 != 12 ? "nd" : (n % 10 == 3 && n % 100 != 13 ? "rd" : "th"))) << " prime is " << prime << std::endl;
    }
    else if (action == "primerange")
    {
        if (argc != 4) {
            std::cerr << "Usage: " << argv[0] << " primerange <start> <end>" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1]) || !primetools::is_numeric(positionals[2])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const mpz_class start(positionals[1].data());
        const mpz_class end(positionals[2].data());

        if (end < start) {
            std::cerr << "Error: End must be greater than or equal to start." << std::endl;
            return 1;
        }

        auto primes = primetools::GetPrimesInRange(start, end);
        std::cout << "Primes in range [" << primetools::TruncateNumber(start) << ", " << primetools::TruncateNumber(end) << "]:" << std::endl;
        for (const auto& prime : primes) {
            std::cout << prime << std::endl;
        }
    }
    else
    {
        std::cerr << "Unknown action: " << action << std::endl;
    }

    return 0;
}