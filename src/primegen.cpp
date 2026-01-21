// primegen
// Precomputes the gaps between primes and outputs them to a file.
// The gaps are encoded as VLE-encoded bytes.
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string_view>

#include <gmpxx.h>

#include "primegenerator.hpp"

std::string
HumanReadableSize(
    const size_t Bytes
) {
    const char* suffixes[] = {"B", "KB", "MB", "GB", "TB"};
    size_t s = 0;
    double count = static_cast<double>(Bytes);
    while (count >= 1024 && s < 4) {
        s++;
        count /= 1024;
    }
    char buffer[50];
    snprintf(buffer, sizeof(buffer), "%.2f %s", count, suffixes[s]);
    return std::string(buffer);
}

int
GeneratePrimeGaps(
    const mpz_class& MaxPrime,
    const bool UseCount,
    const std::string_view OutputFile
) {
    mpz_class last = 0;
    size_t count = 0;
    size_t filesize = 0;
    primetools::PrimeGenerator<mpz_class> generator(9699690);

    std::ofstream ofs(OutputFile.data(), std::ios::binary);
    if (!ofs) {
        std::cerr << "Error: Could not open output file." << std::endl;
        return 1;
    }
    std::cerr << "Generating primes..." << std::endl;
    for (;;) {
        const mpz_class& prime = generator.Next();
        if (prime > MaxPrime) {
            break;
        }
        const mpz_class gap = prime - last;
        last = prime;
        // VLE encode the gap
        mpz_class value = gap;
        while (value > 0) {
            uint8_t byte = value.get_ui() & 0x7F;
            value >>= 7;
            if (value > 0) {
                byte |= 0x80; // Set continuation bit
            }
            ofs.put(static_cast<char>(byte));
            filesize++;
        }
        count++;
        if (UseCount && count >= MaxPrime.get_ui()) {
            break;
        }
        if (count % 100000 == 0) {
            std::cerr << "\rGenerated " << count << " primes, last prime: " << prime << ", file size: " << HumanReadableSize(filesize) << std::flush;
        }
    }
    return 0;
}

int GeneratePrimes(
    const mpz_class& MaxPrime,
    const bool UseCount,
    const std::string_view OutputFile
) {
    size_t count = 0;
    size_t filesize = 0;
    primetools::PrimeGenerator<mpz_class> generator(9699690);

    std::ofstream ofs(OutputFile.data(), std::ios::binary);
    if (!ofs) {
        std::cerr << "Error: Could not open output file." << std::endl;
        return 1;
    }
    std::cerr << "Generating primes..." << std::endl;
    for (;;) {
        const mpz_class& prime = generator.Next();
        if (prime > MaxPrime) {
            break;
        }
        // Output the prime as binary
        uint64_t prime_ui = prime.get_ui();
        ofs.write(reinterpret_cast<const char*>(&prime_ui), sizeof(prime_ui));
        filesize += sizeof(prime_ui);
        count++;
        if (UseCount && count >= MaxPrime.get_ui()) {
            break;
        }
        if (count % 100000 == 0) {
            std::cerr << "\rGenerated " << count << " primes, last prime: " << prime << ", file size: " << HumanReadableSize(filesize) << std::flush;
        }
    }
    return 0;
}

int main(
    int argc,
    char* argv[]
) {
    
    if (argc < 3) {
        std::cerr << "Usage: primegen [options] <output_file>" << std::endl;
        std::cerr << "Options:" << std::endl;
        std::cerr << "  -2 <N>    Generate primes up to 2^N" << std::endl;
        std::cerr << "  -n <N>    Generate primes up to N" << std::endl;
        std::cerr << "  -c <N>    Generate first N primes" << std::endl;
        return 1;
    }

    std::string_view output_file;
    mpz_class max_prime = 0;
    bool use_count = false;
    bool gaps = false;

    for (int i = 1; i < argc; ++i) {
        std::string_view arg = argv[i];
        if (arg == "-2" && i + 1 < argc) {
            mpz_class exp(argv[++i]);
            mpz_class two = 2;
            mpz_pow_ui(max_prime.get_mpz_t(), two.get_mpz_t(), exp.get_ui());
        } else if (arg == "-n" && i + 1 < argc) {
            max_prime = mpz_class(argv[++i]);
        } else if (arg == "-c" && i + 1 < argc) {
            max_prime = mpz_class(argv[++i]);
            use_count = true;
        } else if (arg == "--gaps" || arg == "-g") {
            gaps = true;
        } else {
            output_file = argv[i];
            break;
        }
    }

    if (output_file.empty() || max_prime == 0) {
        std::cerr << "Error: Invalid arguments." << std::endl;
        return 1;
    }    

    int result = gaps ?
        GeneratePrimeGaps(max_prime, use_count, output_file) :
        GeneratePrimes(max_prime, use_count, output_file);
    if (result != 0) {
        return result;
    }

    std::cerr << std::endl << "Finished generating primes." << std::endl;
    std::cerr << "Output file size: " << HumanReadableSize(std::filesystem::file_size(output_file)) << std::endl;
    return 0;
}