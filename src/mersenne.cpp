
#include <iostream>
#include <string>
#include <string_view>

#include "gmpxx.h"

#include "factordb.hpp"
#include "factorise.hpp"
#include "util.hpp"

static mpz_class gCurrentPower = 0;
static mpz_class gCurrentMersenne = 0;
static size_t gLastLogLength = 0;

void
Log(
    std::string_view Message
)
{
    std::cout << '\r' << std::string(gLastLogLength, ' ') << std::flush;
    std::cout << '\r' << "2" << primetools::ToSuperScript(gCurrentPower) << "-1 = " << gCurrentMersenne << " = " << Message << std::flush;
    gLastLogLength = 3 + primetools::ToString(gCurrentPower).size() + 7 + primetools::ToString(gCurrentMersenne).size() + 3 + Message.size();
}

int main(
    int argc,
    char* argv[]
) {
    if (argc < 2) {
        std::cerr << "Usage: mersenne [options] <max_power>" << std::endl;
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
            std::cerr << "Usage: mersenne [options] <max_power>" << std::endl;
            std::cerr << "Options:" << std::endl;
            std::cerr << "  -d <path>   Path to factor database" << std::endl;
            std::cerr << "  -h          Show this help message" << std::endl;
            return 0;
        } else if (max_power == 0) {
            max_power = argv[i];
        }
    }

    primetools::PreCacheWheelAndPrimesInNewThread();

    std::cerr << "Finding prime factors for Mersenne numbers up to 2" << primetools::ToSuperScript(max_power) << "-1" << std::endl;
    for (mpz_class p = 2; p <= max_power; p++) {
        mpz_class mersenne;
        mpz_ui_pow_ui(mersenne.get_mpz_t(), 2, p.get_ui());
        mersenne -= 1;

        gCurrentPower = p;
        gCurrentMersenne = mersenne;

        std::cout << "2" << primetools::ToSuperScript(p) << "-1 = " << mersenne << " = ";

        // First check if it is prime
        if (primetools::IsPrime(mersenne)) {
            std::cout << mersenne << std::endl;
            continue;
        }

        const auto factors = primetools::Factorise(mersenne, 0, db_path, Log);
        if (factors) {
            Log(factors->GetString());   
        }
        std::cout << std::endl;
    }
}