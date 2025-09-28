
#include <stdint.h>
#include <gmpxx.h>

#include <limits>
#include <iostream>
#include <string>

#include "factorise.hpp"
#include "util.hpp"

void
FindClosePrimes(
    void
)
{
    mpz_class a, b;
    mpz_class end;
    mpz_class smallest, smallestdiff = 0;

    // Set a to 2^1023
    mpz_ui_pow_ui(a.get_mpz_t(), 2, 1023);

    // Set end to 2^1024
    mpz_ui_pow_ui(end.get_mpz_t(), 2, 1024);

    // Set a to the first prime greater than 2^1023
    mpz_nextprime(a.get_mpz_t(), a.get_mpz_t());

    // Releatedly look for primes with the smallest distance
    // between them
    while (a < end) {
        // Find the next prime
        mpz_nextprime(b.get_mpz_t(), a.get_mpz_t());

        // Check if the difference is smaller than the smallest
        if (smallestdiff == 0 || (b - a) < smallestdiff) {
            smallestdiff = b - a;
            smallest = a;
            std::cout << "Found a smaller pair of primes: "
                      << smallest << ", " << b <<
                      " (" << smallestdiff << ")" << std::endl;
        }
        a = b;
    }
}

const std::string
TruncateNumber(
    const mpz_class& Number,
    const size_t StartDigits = 5,
    const size_t EndDigits = 3
)
{
    std::string str = Number.get_str();
    if (str.length() <= StartDigits + EndDigits) {
        return str;
    }
    return str.substr(0, StartDigits) + "..." + str.substr(str.length() - EndDigits);
}

void
OutputFactors(
    const std::optional<std::pair<mpz_class, mpz_class>>& Factors
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

    std::string action = argv[1];

    if (action == "factorise" || action == "factorize")
    {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " factorise <number>" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);

        auto result = primetools::Factorise(n);
        OutputFactors(result);
    }
    else if (action == "findcloseprimes")
    {
        FindClosePrimes();
    }
    else if (action == "isprime")
    {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " isprime <number>" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);

        bool is_prime = primetools::isprime(n);
        std::cout << TruncateNumber(n) << (is_prime ? " is prime." : " is not prime.") << std::endl;
    }
    else if (action == "fermat")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " fermat <number>" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);
        const size_t max_iterations = argc == 4 ? std::stoul(argv[3]) : std::numeric_limits<size_t>::max();

        std::cout << "Fermat factorization of " << TruncateNumber(n) << std::endl;

        auto result = primetools::FermatFactorisation(n, 0, max_iterations);
        OutputFactors(result);
    }
    else if (action == "fermat2" || action == "ffa2")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " fermat <number>" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);
        const size_t max_iterations = argc == 4 ? std::stoul(argv[3]) : std::numeric_limits<size_t>::max();

        std::cout << "Fermat factorization of " << TruncateNumber(n) << " using FFA-2" << std::endl;

        auto result = primetools::FermatFactorisationAlgorithm2(n, max_iterations);
        OutputFactors(result);
    }
    else if (action == "fermat4" || action == "mffv4")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " fermat <number>" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);
        const size_t max_iterations = argc == 4 ? std::stoul(argv[3]) : std::numeric_limits<size_t>::max();

        std::cout << "Fermat factorization of " << TruncateNumber(n) << " using MFFV4" << std::endl;

        auto result = primetools::ModifiedFermatFactorisation4(n, max_iterations);
        OutputFactors(result);
    }
    else if (action == "fmmod20pc")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " fmmod20pc <number>" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);
        const size_t max_iterations = argc == 4 ? std::stoul(argv[3]) : std::numeric_limits<size_t>::max();

        std::cout << "Fermat factorization of " << TruncateNumber(n) << " using FMMod20Precomp" << std::endl;

        auto result = primetools::FMMod20Precomp(n, max_iterations);
        OutputFactors(result);
    }
    else if (action == "pollardsrho")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " pollardsrho <number>" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);

        std::cout << "Pollard's rho factorization of " << TruncateNumber(n) << std::endl;

        auto result = primetools::PollardsRho(n);
        OutputFactors(result);
    }
    else if (action == "brentpollardsrho")
    {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " brentpollardsrho <number>" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);

        auto result = primetools::BrentPollardsRho(n);
        OutputFactors(result);
    }
    else if (action == "pollardp1")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " pollardp1 <number> <bound>" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);
        const size_t bound = argc == 4 ? std::stoul(argv[3]) : (size_t)1 << 20;

        auto result = primetools::PollardsPMinus1(n, bound);
        OutputFactors(result);
    }
    else if (action == "squfof")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " squfof <number>" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);

        auto result = primetools::SQUFOF(n);
        OutputFactors(result);
    }
    else if (action == "random" || action == "fuckit")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " randomprimefactorization <number> [max_iterations]" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);

        std::cout << "Random prime factorization of " << TruncateNumber(n) << std::endl;
        
        const size_t base = argc >= 4 ? std::stoul(argv[3]) : 1;
        const size_t max_iterations = argc == 5 ? std::stoul(argv[4]) : std::numeric_limits<size_t>::max();

        auto result = primetools::TrialDivisionRandom(n, base, max_iterations);
        OutputFactors(result);
    }
    else if (action == "wheel30")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " wheel30 <number> [max_iterations]" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);
        const size_t max_iterations = argc == 4 ? std::stoul(argv[3]) : std::numeric_limits<size_t>::max();

        auto result = primetools::TrialDivisionWheel30(n, max_iterations);
        OutputFactors(result);
    }
    else if (action == "wheel210")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " wheel30 <number> [max_iterations]" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);
        const size_t max_iterations = argc == 4 ? std::stoul(argv[3]) : std::numeric_limits<size_t>::max();

        auto result = primetools::TrialDivisionWheel210(n, max_iterations);
        OutputFactors(result);
    }
    else if (action == "trialdivision" || action == "linear" || action == "trial" || action == "td")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " trialdivision <number> [base] [max_iterations]" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);
        const size_t base = argc >= 4 ? std::stoul(argv[3]) : 1;
        const size_t max_iterations = argc == 5 ? std::stoul(argv[4]) : std::numeric_limits<size_t>::max();

        auto result = primetools::TrialDivisionLinear(n, base, max_iterations);
        OutputFactors(result);
    }
    else if (action == "bitflip")
    {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " trialdivision <number> [max_iterations]" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);
        const size_t max_iterations = argc == 5 ? std::stoul(argv[4]) : std::numeric_limits<size_t>::max();

        auto result = primetools::TrialDivisionBitflip(n, max_iterations);
        OutputFactors(result);
    }
    else if (action == "findcloseprimes")
    {
        FindClosePrimes();
    }
    else if (action == "calculatefermatiterations")
    {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " calculatefermatiterations <number>" << std::endl;
            return 1;
        }

        const mpz_class n(argv[2]);
        size_t iterations = primetools::CalculateFermatIterations(n);
        std::cout << "Fermat iterations for " << TruncateNumber(n) << ": " << iterations << std::endl;
    }
    else
    {
        std::cerr << "Unknown action: " << action << std::endl;
    }

    return 0;
}