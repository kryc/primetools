#ifndef PRIMETOOLS_FACTORDB_HPP
#define PRIMETOOLS_FACTORDB_HPP

#include <algorithm>
#include <cinttypes>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <optional>
#include <span>
#include <sys/mman.h>
#include <string_view>
#include <vector>

#include "factors.hpp"
#include "util.hpp"

namespace primetools {

template<typename T>
class FactorDB {
public:
    FactorDB(
        const std::string_view path = ""
    ) : m_CachePath(path) {
        if (m_CachePath.empty()) {
            return;
        }

        if (!std::filesystem::exists(m_CachePath)) {
            std::filesystem::create_directories(m_CachePath);
        }
    };

    ~FactorDB() {
        // Close();
    }

    const bool
    IsOpen(
        void
    ) const {
        return !m_CachePath.empty();
    }

    const std::filesystem::path
    GetPath(
        void
    ) const {
        return m_CachePath;
    }

    const bool
    Clear(
        void
    ) {
        try {
            std::filesystem::remove_all(m_CachePath);
            std::filesystem::create_directories(m_CachePath);
            return true;
        } catch (...) {
            return false;
        }
    }

    const uint32_t
    GetFactorKey(
        const T Value
    ) const {
        // The key is the lowest five bytes of the number
        if constexpr (std::is_same_v<T, mpz_class>) {
            auto uival = Value.get_ui();
            return static_cast<uint32_t>(uival & 0xFFFFFFFFFF);
        } else {
            return static_cast<uint32_t>(Value & 0xFFFFFFFFFF);
        }
    }

    const std::filesystem::path
    GetFactorPathBase(
        const T Value,
        const std::string Suffix
    ) const {
        const uint32_t key = GetFactorKey(Value);
        const std::string key_hex = primetools::GetHex(key, 6);
        return m_CachePath / ("factors_" + key_hex + Suffix);
    }

    const std::filesystem::path
    GetFactorIndex(
        const T Value
    ) const {
        return GetFactorPathBase(Value, ".idx");
    }

    const std::filesystem::path
    GetFactorData(
        const T Value
    ) const {
        return GetFactorPathBase(Value, ".dat");
    }

    void AddFactors(
        const PrimeFactors<T> Factors
    ) {
        AddFactorsInternal(Factors);
    }

    std::optional<PrimeFactors<T>>
    GetFactors(
        const T& Product
    ) const {
        auto indexPath = GetFactorIndex(Product);
        auto dataPath = GetFactorData(Product);
        if (!std::filesystem::exists(indexPath) || !std::filesystem::exists(dataPath)) {
            return std::nullopt;
        }

        // Open the index file and read offsets
        std::ifstream indexFile(indexPath, std::ios::binary);
        if (!indexFile.is_open()) {
            throw std::runtime_error("Failed to open FactorDB index file for reading: " + indexPath.string());
        }
        const size_t indexFileSize = std::filesystem::file_size(indexPath);
        std::vector<size_t> offsets;
        offsets.resize(indexFileSize / sizeof(size_t));
        indexFile.read(reinterpret_cast<char*>(offsets.data()), indexFileSize);
        indexFile.close();

        // Open the data file and read factors
        std::ifstream data_file(dataPath, std::ios::binary);
        if (!data_file.is_open()) {
            throw std::runtime_error("Failed to open FactorDB data file for reading: " + dataPath.string());
        }

        PrimeFactors<T> factors;
        T remainder;
        for (const size_t off : offsets) {
            data_file.seekg(off);
            // Read product as VLE, byte by byte from the data file
            T nextProduct = 0;
            uint8_t byte = 0;
            size_t shift = 0;
            do {
                data_file.read(reinterpret_cast<char*>(&byte), 1);
                nextProduct |= static_cast<T>(byte & 0x7F) << shift;
                shift += 7;
            } while (byte & 0x80);
            // Read number of factors as VLE. We will always need this
            size_t numFactors = 0;
            shift = 0;
            do {
                data_file.read(reinterpret_cast<char*>(&byte), 1);
                numFactors |= static_cast<size_t>(byte & 0x7F) << shift;
                shift += 7;
            } while (byte & 0x80);
            if (nextProduct != Product) {
                continue; // Not the product we are looking for
            }
            // We have found the product, read the factors
            remainder = Product;
            for (size_t i = 0; i < numFactors; ++i) {
                T factor = 0;
                shift = 0;
                do {
                    data_file.read(reinterpret_cast<char*>(&byte), 1);
                    factor |= static_cast<T>(byte & 0x7F) << shift;
                    shift += 7;
                } while (byte & 0x80);
                while (remainder % factor == 0) {
                    factors.AddFactor(factor);
                    remainder /= factor;
                }
            }
            return factors;
        }
        data_file.close();
        return factors;
    }

    const size_t
    GetCount(
        void
    ) const {
        size_t total = 0;
        for (const auto& entry : std::filesystem::directory_iterator(m_CachePath)) {
            if (entry.path().extension() == ".idx") {
                const size_t file_size = std::filesystem::file_size(entry.path());
                total += file_size / sizeof(size_t);
            }
        }
        return total;
    }

    void PrintStats(
        void
    ) const {
        
    }
private:
    void AddFactorsInternal(
        const PrimeFactors<T> Factors
    ) {
        T product = Factors.Product();
        auto indexPath = GetFactorIndex(product);
        auto dataPath = GetFactorData(product);
        // Save the current size of the data file
        size_t currentOffset = std::filesystem::exists(dataPath)
            ? std::filesystem::file_size(dataPath)
            : 0;
        // Get the serialized factor data
        std::vector<uint8_t> serializedData = Factors.Serialize();
        // Open the data file for appending
        std::ofstream data_file(dataPath, std::ios::binary | std::ios::app);
        if (!data_file.is_open()) {
            throw std::runtime_error("Failed to open FactorDB data file for writing: " + dataPath.string());
        }
        data_file.write(reinterpret_cast<const char*>(serializedData.data()), serializedData.size());
        data_file.close();
        // Now update the index file with the offset of the new data
        std::ofstream indexFile(indexPath, std::ios::binary | std::ios::app);
        if (!indexFile.is_open()) {
            throw std::runtime_error("Failed to open FactorDB index file for writing: " + indexPath.string());
        }
        indexFile.write(reinterpret_cast<const char*>(&currentOffset), sizeof(currentOffset));
        indexFile.close();
    }

    std::filesystem::path m_CachePath;
};

} // namespace primetools

#endif // PRIMETOOLS_FACTORDB_HPP