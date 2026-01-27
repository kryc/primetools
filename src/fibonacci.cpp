
#include <iostream>
#include <string>
#include <string_view>

#include "gmpxx.h"

#include "factordb.hpp"
#include "factorise.hpp"
#include "util.hpp"

static mpz_class gCurrentFibonacci = 0;
static size_t gLastLogLength = 0;

void
Log(
    std::string_view Message
)
{
    std::cout << '\r' << std::string(gLastLogLength, ' ') << std::flush;
    std::cout << '\r' << gCurrentFibonacci << " = " << Message << std::flush;
    gLastLogLength = primetools::ToString(gCurrentFibonacci).size() + 3 + Message.size();
}

int main(
    int argc,
    char* argv[]
) {
    if (argc < 2) {
        std::cerr << "Usage: fibonacci [options] <max_power>" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "  -d <path>   Path to factor database" << std::endl;
        std::cerr << "  -h          Show this help message" << std::endl;
        return 1;
    }

    std::string_view db_path;
    mpz_class max_power;

    for (size_t i = 1; i < argc; i++) {
        std::string_view arg = argv[i];
        if ((arg == "-d" || arg == "--db") && i + 1 < argc) {
            db_path = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            std::cerr << "Usage: fibonacci [options] <max_power>" << std::endl;
            std::cerr << "Options:" << std::endl;
            std::cerr << "  -d <path>   Path to factor database" << std::endl;
            std::cerr << "  -h          Show this help message" << std::endl;
            return 0;
        } else if (max_power == 0) {
            max_power = argv[i];
        }
    }

    std::cerr << "Finding prime factors for Fibonacci numbers up to 2" << primetools::ToSuperScript(max_power) << std::endl;
    mpz_class end;
    mpz_ui_pow_ui(end.get_mpz_t(), 2, max_power.get_ui());
    mpz_class last = 0;
    gCurrentFibonacci = 1;
    while (gCurrentFibonacci <= end) {
        gCurrentFibonacci = last + gCurrentFibonacci;
        last = gCurrentFibonacci - last;

        std::cout << gCurrentFibonacci << " = ";

        // First check if it is prime
        if (primetools::IsPrime(gCurrentFibonacci)) {
            std::cout << gCurrentFibonacci << std::endl;
            continue;
        }

        const auto factors = primetools::Factorise(gCurrentFibonacci, 0, db_path, Log);
        if (factors) {
            Log(factors->GetString());   
        }
        std::cout << std::endl;
    }
}