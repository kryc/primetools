//
// FactorDB
// A simple to manage the factor database
//

#include <iostream>
#include <iomanip>
#include <optional>
#include <span>
#include <string>
#include <string_view>
#include <thread>

#include <gmpxx.h>

#include "factordb.hpp"
#include "factors.hpp"

using namespace primetools;

static const std::string_view HELP_STRING = R"(
Usage: factordb [path] [action] [options] args
Actions:
    get <number>        Retrieve factors for <number>
    add <factors>       Add factors to the database
    count               Get the number of products in the database
)";

int main(
    int argc,
    char* argv[]
) {

    if (argc < 3) {
        std::cerr << HELP_STRING << std::endl;
        return 1;
    }

    std::string_view db_path;
    std::string_view action;
    std::vector<std::string_view> positionals;
    std::map<std::string_view, std::string_view> flags;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (db_path.empty()) {
            db_path = arg;
        } else if (action.empty()) {
            action = arg;
        } else if (arg == "--help" || arg == "-h") {
            std::cout << HELP_STRING << std::endl;
            return 0;
        } else if (arg.rfind("--", 0) == 0) {
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

    FactorDB<mpz_class> db(db_path);
    if (action == "get") {
        if (positionals.size() != 1) {
            std::cerr << "Usage: factordb <path> get <number>" << std::endl;
            return 1;
        }

        if (!primetools::is_numeric(positionals[0])) {
            std::cerr << "Error: Input is not a valid number." << std::endl;
            return 1;
        }

        const mpz_class n(positionals[0].data());
        auto factors_opt = db.GetFactors(n);
        if (factors_opt) {
            PrimeFactors<mpz_class> factors = factors_opt.value();
            std::cout << "Factors of " << n << ": " << factors.GetString() << std::endl;
        } else {
            std::cout << "No factors found for " << n << " in the database." << std::endl;
        }
    } else if (action == "add") {
        if (positionals.size() < 2) {
            std::cerr << "Usage: factordb <path> add <factor1> <factor2> ..." << std::endl;
            return 1;
        }

        PrimeFactors<mpz_class> factors;
        for (const auto& factor_str : positionals) {
            if (!primetools::is_numeric(factor_str)) {
                std::cerr << "Error: Factor '" << factor_str << "' is not a valid number." << std::endl;
                return 1;
            }
            const mpz_class factor(factor_str.data());
            if (factor < 2) {
                std::cerr << "Error: Factor '" << factor_str << "' must be greater than 1." << std::endl;
                return 1;
            }
            if (!primetools::isprime(factor)) {
                std::cerr << "Error: Factor '" << factor_str << "' is not prime." << std::endl;
                return 1;
            }
            factors.AddFactor(factor);
        }
        if (factors.Size() == 0) {
            std::cerr << "Error: No valid factors provided." << std::endl;
            return 1;
        } else if (factors.Size() == 1) {
            std::cerr << "Error: At least two factors are required to form a product." << std::endl;
            return 1;
        }
        // Check if it is already in the database
        if (db.GetFactors(factors.Product())) {
            std::cerr << "Error: Factors for product " << factors.Product() << " already exist in the database." << std::endl;
            return 1;
        }
        db.AddFactors(factors);
        std::cout << "Added factors to the database: " << factors.GetString() << std::endl;
    } else if (action == "count") {
        std::cout << db.GetCount() << std::endl;
        return 1;
    } else if (action == "rebuildindex") {
        size_t new_index_bits = 32;
        if (flags.find("bits") != flags.end()) {
            new_index_bits = std::stoul(std::string(flags["bits"]));
            if (new_index_bits != 32 && new_index_bits != 64) {
                std::cerr << "Error: Index bits must be either 32 or 64." << std::endl;
                return 1;
            }
        }
        db.RebuildIndex(new_index_bits);
        std::cout << "Rebuilt index files to " << new_index_bits << "-bit." << std::endl;
    } else {
        std::cerr << "Unknown action: " << action << std::endl;
        return 1;
    }
}