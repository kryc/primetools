
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

        mpz_class n(argv[2]);

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

        mpz_class n(argv[2]);

        bool is_prime = primetools::isprime(n);
        std::cout << n << (is_prime ? " is prime." : " is not prime.") << std::endl;
    }
    else if (action == "pollardsrho")
    {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " pollardsrho <number>" << std::endl;
            return 1;
        }

        mpz_class n(argv[2]);

        auto result = primetools::PollardsRho(n);
        OutputFactors(result);
    }
    else
    {
        std::cerr << "Unknown action: " << action << std::endl;
    }

    return 0;
}