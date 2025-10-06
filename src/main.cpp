#include <iostream>
#include <limits>
#include <map>
#include <string>

#include <stdint.h>
#include <gmpxx.h>

#include "factorise.hpp"
#include "util.hpp"

template <typename T>
void
OutputFactors(
    const std::optional<std::pair<T, T>>& Factors
)
{
    if (Factors) {
        std::cout << Factors->first << ", " << Factors->second << std::endl;
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

    std::map<std::string, std::string> flags;
    std::vector<std::string> positionals;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg.rfind("--", 0) == 0) {
            size_t eq_pos = arg.find('=');
            if (eq_pos != std::string::npos) {
                std::string key = arg.substr(2, eq_pos - 2);
                std::string value = arg.substr(eq_pos + 1);
                flags[key] = value;
            } else {
                std::string key = arg.substr(2);
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

    const std::string action = positionals[0];

    if (action == "factorise" || action == "factorize")
    {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " factorise <number>" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const mpz_class n(positionals[1]);

        auto result = primetools::Factorise(n);
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

        const mpz_class n(positionals[1]);
        const size_t max = flags.find("max-iterations") != flags.end() ? std::stoul(flags["max-iterations"]) : std::numeric_limits<size_t>::max();
        const bool fermat = flags.find("fermat") != flags.end() || flags.find("ffa") != flags.end();
        const bool mffv2 = flags.find("mffv2") != flags.end();
        const bool mffv4 = flags.find("mffv4") != flags.end();

        if (fermat) {
            std::cout << "Using standard Fermat factorization of " << primetools::TruncateNumber(n) << std::endl;
            auto result = primetools::FermatFactorisation(n, 0, max);
            OutputFactors(result);
        }
        else if (mffv2) {
            std::cout << "Using Fermat Factorization Algorithm 2 of " << primetools::TruncateNumber(n) << std::endl;
            auto result = primetools::FermatFactorisationAlgorithm2(n, max);
            OutputFactors(result);
        }
        else if (mffv4) {
            std::cout << "Using Modified Fermat Factorization V4 of " << primetools::TruncateNumber(n) << std::endl;
            auto result = primetools::ModifiedFermatFactorisation4(n, max);
            OutputFactors(result);
        }
        else {
            std::cout << "Using FMMod20Precomp of " << primetools::TruncateNumber(n) << std::endl;
            auto result = primetools::FMMod20Precomp(n, max);
            OutputFactors(result);
        }
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

        const mpz_class n(positionals[1]);
        const bool pollard = flags.find("pollard") != flags.end();

        if (pollard) {
            std::cout << "Pollard's rho factorization of " << primetools::TruncateNumber(n) << std::endl;
            auto result = primetools::PollardsRho(n);
            OutputFactors(result);
        }
        else {
            std::cout << "Brent-Pollard's rho factorization of " << primetools::TruncateNumber(n) << std::endl;
            auto result = primetools::BrentPollardsRho(n);
            OutputFactors(result);
        }
    }
    else if (action == "pminus1" || action == "p-1" || action == "pollardp1" || action == "p1" || action == "pollardp-1")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " pminus1 <number> <bound>" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const mpz_class n(positionals[1]);
        const size_t bound = positionals.size() >= 3 ? std::stoul(positionals[2]) : (1 << 20);
        auto result = primetools::PollardsPMinus1(n, bound);
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

        const mpz_class n(positionals[1]);
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
            modulus = flags.find("modulus") != flags.end() ? std::stoul(flags["modulus"]) : std::stoul(flags["wheel"]);
        }

        const mpz_class n(positionals[1]);
        const mpz_class start = flags.find("start") != flags.end() ? mpz_class(flags["start"]) : mpz_class(1);
        const mpz_class end = flags.find("end") != flags.end() ? mpz_class(flags["end"]) : mpz_class(0);
        const bool noguess = flags.find("no-guess") != flags.end();
        const size_t bits = flags.find("bits") != flags.end() ? std::stoul(flags["bits"]) : 0;
        const size_t threads = flags.find("threads") != flags.end() ? std::stoul(flags["threads"]) : 1;
        const size_t blocksize = flags.find("blocksize") != flags.end() ? std::stoul(flags["blocksize"]) : 0;
        const uint64_t seed = flags.find("seed") != flags.end() ? std::stoull(flags["seed"]) : 0;
        const size_t max_iterations = flags.find("max-iterations") != flags.end() ? std::stoul(flags["max-iterations"]) : std::numeric_limits<size_t>::max();

        std::cout << "Random prime factorization of " << primetools::TruncateNumber(n) << std::endl;

        if (threads == 0 || threads > 1) {
            auto result = primetools::TrialDivisionRandomMT(n, threads, blocksize, !noguess, bits, start, end, modulus);
            OutputFactors(result);
        }
        else {
            auto result = primetools::TrialDivisionRandom(n, !noguess, bits, start, end, seed, modulus, max_iterations);
            OutputFactors(result);
        }
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
            modulus = flags.find("modulus") != flags.end() ? std::stoul(flags["modulus"]) : std::stoul(flags["wheel"]);
        }

        const mpz_class n(positionals[1]);
        const mpz_class start = flags.find("start") != flags.end() ? mpz_class(flags["start"]) : mpz_class(0);
        const mpz_class end = flags.find("end") != flags.end() ? mpz_class(flags["end"]) : mpz_class(0);
        const bool noguess = flags.find("no-guess") != flags.end();
        const size_t bits = flags.find("bits") != flags.end() ? std::stoul(flags["bits"]) : 0;
        const bool simd = flags.find("simd") != flags.end() || flags.find("use-simd") != flags.end();
        const size_t threads = flags.find("threads") != flags.end() ? std::stoul(flags["threads"]) : 1;
        const size_t blocksize = flags.find("blocksize") != flags.end() ? std::stoul(flags["blocksize"]) : 0;

        if (simd) {
            const size_t max_iterations = flags.find("max-iterations") != flags.end() ? std::stoul(flags["max-iterations"]) : std::numeric_limits<size_t>::max();
            auto result = primetools::TrialDivisionSimd(n, max_iterations);
            OutputFactors(result);
        }
        else if (threads == 0 || threads > 1) {
            auto result = primetools::TrialDivisionMT(n, threads, blocksize, !noguess, bits, start, end, modulus);
            OutputFactors(result);
        }
        else {
            if (n.fits_ulong_p()) {
                auto result = primetools::TrialDivision<uint64_t>(n.get_ui(), modulus, !noguess, bits, start.get_ui(), end.get_ui());
                OutputFactors(result);
            }
            else
            {
                auto result = primetools::TrialDivision<mpz_class>(n, modulus, !noguess, bits, start, end);
                OutputFactors(result);
            }  
        }
    }
    else if (action == "bitflip")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " trialdivision <number> [max_iterations]" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[1])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const mpz_class n(positionals[1]);
        const size_t max_iterations = argc == 5 ? std::stoul(argv[4]) : std::numeric_limits<size_t>::max();

        auto result = primetools::TrialDivisionBitflip(n, max_iterations);
        OutputFactors(result);
    }
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

        const mpz_class n(positionals[1]);

        bool is_prime = primetools::isprime(n);
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

        const mpz_class n(positionals[1]);
        size_t iterations = primetools::CalculateFermatIterations(n);
        std::cout << "Fermat iterations for " << primetools::TruncateNumber(n) << ": " << iterations << std::endl;
    }
    else
    {
        std::cerr << "Unknown action: " << action << std::endl;
    }

    return 0;
}