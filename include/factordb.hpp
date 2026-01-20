#ifndef PRIMETOOLS_FACTORDB_HPP
#define PRIMETOOLS_FACTORDB_HPP

#include <algorithm>
#include <cinttypes>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
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
    constexpr static std::string_view kIndexFile64Ext = ".idx";
    constexpr static std::string_view kIndexFile32Ext = ".i32";
    constexpr static std::string_view kDataFileExt = ".dat";
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
            return static_cast<uint32_t>(uival & 0xFFFFF);
        } else {
            return static_cast<uint32_t>(Value & 0xFFFFF);
        }
    }

    const std::filesystem::path
    GetFactorPathBase(
        const T Value,
        const std::string_view Suffix
    ) const {
        const uint32_t key = GetFactorKey(Value);
        const std::string key_hex = primetools::GetHex(key, 5);
        return m_CachePath / ("factors_" + key_hex + std::string(Suffix));
    }

    const std::filesystem::path
    GetFactorIndex(
        const T Value
    ) const {
        // Check if 32 or 64 bit index exists
        auto path64 = GetFactorPathBase(Value, kIndexFile64Ext);
        if (std::filesystem::exists(path64)) {
            return path64;
        }
        return GetFactorPathBase(Value, kIndexFile32Ext);
    }

    static const size_t
    GetIndexBitsFromPath(
        const std::filesystem::path& Path
    ) {
        if (Path.extension() == kIndexFile32Ext) {
            return 32;
        } else {
            return 64;
        }
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
        // Check if already exists
        if (GetFactors(Factors.Product())) {
            return;
        }
        AddFactorsInternal(Factors);
    }

    template <typename U>
    std::optional<PrimeFactors<T>>
    FindFactor(
        const T& N,
        std::ifstream& DataFile,
        const size_t DataFileSize,
        std::ifstream& IndexFile,
        const size_t IndexFileSize
    ) const {
        (void)DataFileSize;

        // Index file consists of U (uint32/uint64_t) offsets into the data file.
        // Requires the index to be sorted by product (see SortIndexInternal).
        const size_t numEntries = IndexFileSize / sizeof(U);
        if (numEntries == 0) {
            return std::nullopt;
        }

        auto readVleT = [&](std::ifstream& file) {
            T value = 0;
            uint8_t byte = 0;
            size_t shift = 0;
            do {
                file.read(reinterpret_cast<char*>(&byte), 1);
                value |= static_cast<T>(byte & 0x7F) << shift;
                shift += 7;
            } while (byte & 0x80);
            return value;
        };

        auto readVleSizeT = [&](std::ifstream& file) {
            size_t value = 0;
            uint8_t byte = 0;
            size_t shift = 0;
            do {
                file.read(reinterpret_cast<char*>(&byte), 1);
                value |= static_cast<size_t>(byte & 0x7F) << shift;
                shift += 7;
            } while (byte & 0x80);
            return value;
        };

        auto readOffsetAtIndex = [&](const size_t index) {
            U offset;
            IndexFile.clear();
            IndexFile.seekg(static_cast<std::streamoff>(index * sizeof(U)));
            IndexFile.read(reinterpret_cast<char*>(&offset), sizeof(U));
            return offset;
        };

        auto readProductAtOffset = [&](const U offset) {
            DataFile.clear();
            DataFile.seekg(static_cast<std::streamoff>(offset));
            return readVleT(DataFile);
        };

        // lower_bound search for N.
        size_t left = 0;
        size_t right = numEntries;
        while (left < right) {
            const size_t mid = left + (right - left) / 2;
            const U offset = readOffsetAtIndex(mid);
            const T product = readProductAtOffset(offset);
            if (product < N) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }

        if (left >= numEntries) {
            return std::nullopt;
        }

        const U offset = readOffsetAtIndex(left);
        DataFile.clear();
        DataFile.seekg(static_cast<std::streamoff>(offset));
        const T product = readVleT(DataFile);
        if (product != N) {
            return std::nullopt;
        }

        const size_t numFactors = readVleSizeT(DataFile);
        PrimeFactors<T> factors;
        for (size_t j = 0; j < numFactors; ++j) {
            const T factor = readVleT(DataFile);
            factors.AddFactor(factor);
        }
        return factors;
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
        const size_t bits = GetIndexBitsFromPath(indexPath);

        auto dataFile = std::ifstream(dataPath, std::ios::binary);
        if (!dataFile.is_open()) {
            throw std::runtime_error("Failed to open FactorDB data file for reading: " + dataPath.string());
        }
        const size_t dataFileSize = std::filesystem::file_size(dataPath);

        std::optional<PrimeFactors<T>> result;
        if (bits == 32) {
            result = FindFactor<uint32_t>(Product, dataFile, dataFileSize, indexFile, indexFileSize);
        } else {
            result = FindFactor<uint64_t>(Product, dataFile, dataFileSize, indexFile, indexFileSize);
        }
        // Close files
        dataFile.close();
        indexFile.close();
        return result;
    }

    const size_t
    GetCount(
        void
    ) const {
        size_t total = 0;
        for (const auto& entry : std::filesystem::directory_iterator(m_CachePath)) {
            if (entry.path().extension() == ".idx") {
                const size_t file_size = std::filesystem::file_size(entry.path());
                total += file_size / sizeof(uint64_t);
            } else if (entry.path().extension() == ".i32") {
                const size_t file_size = std::filesystem::file_size(entry.path());
                total += file_size / sizeof(uint32_t);
            }
        }
        return total;
    }

    const bool
    RebuildIndex(
        const size_t NewIndexBits = 32
    ) const {
        // Check bits are 32 or 64
        if (NewIndexBits != 32 && NewIndexBits != 64) {
            return false;
        }
        // First, delete all existing index files
        for (const auto& entry : std::filesystem::directory_iterator(m_CachePath)) {
            if (entry.path().extension() == kIndexFile32Ext || entry.path().extension() == kIndexFile64Ext) {
                std::filesystem::remove(entry.path());
            }
        }
        // Now, rebuild the index files from the data files
        for (const auto& entry : std::filesystem::directory_iterator(m_CachePath)) {
            if (entry.path().extension() != ".dat") {
                continue;
            }
            std::ifstream data_file(entry.path(), std::ios::binary);
            if (!data_file.is_open()) {
                throw std::runtime_error("Failed to open FactorDB data file for reading: " + entry.path().string());
            }
            auto index_file_path = std::filesystem::path(entry.path()).replace_extension(GetIndexExtensionForBits(NewIndexBits));
            std::ofstream index_file(
                index_file_path,
                std::ios::binary | std::ios::app
            );
            if (!index_file.is_open()) {
                throw std::runtime_error("Failed to open FactorDB index file for writing: " + index_file_path.string());
            }
            size_t currentOffset = 0;
            while (data_file.peek() != EOF) {
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
                // Skip the factors
                for (size_t i = 0; i < numFactors; ++i) {
                    do {
                        data_file.read(reinterpret_cast<char*>(&byte), 1);
                    } while (byte & 0x80);
                }
                // Write the offset to the index file
                if (NewIndexBits == 32) {
                    uint32_t offset32 = static_cast<uint32_t>(currentOffset);
                    index_file.write(reinterpret_cast<const char*>(&offset32), sizeof(offset32));
                } else {
                    uint64_t offset64 = static_cast<uint64_t>(currentOffset);
                    index_file.write(reinterpret_cast<const char*>(&offset64), sizeof(offset64));
                }
                // Update the current offset
                currentOffset = data_file.tellg();;
            }
            data_file.close();
            index_file.close();
        }
        // Now sort all index files
        for (const auto& entry : std::filesystem::directory_iterator(m_CachePath)) {
            if (entry.path().extension() == kIndexFile32Ext || entry.path().extension() == kIndexFile64Ext) {
                auto data_file_path = std::filesystem::path(entry.path()).replace_extension(kDataFileExt);
                SortIndex(entry.path(), data_file_path);
            }
        }
        // Now we are done
        return true;
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

        // Index entries are offsets into the data file; width must match the index file extension.
        if (indexPath.extension() == kIndexFile32Ext) {
            if (currentOffset > static_cast<size_t>(std::numeric_limits<uint32_t>::max())) {
                throw std::runtime_error(
                    "FactorDB data file too large for 32-bit index (" + indexPath.string() + "). Rebuild index as 64-bit."
                );
            }
            const uint32_t offset32 = static_cast<uint32_t>(currentOffset);
            indexFile.write(reinterpret_cast<const char*>(&offset32), sizeof(offset32));
        } else {
            const uint64_t offset64 = static_cast<uint64_t>(currentOffset);
            indexFile.write(reinterpret_cast<const char*>(&offset64), sizeof(offset64));
        }
        indexFile.close();
        // Sort the index file
        SortIndex(indexPath, dataPath);
    }

    const size_t
    GetIndexFileCount(
        const std::filesystem::path& IndexPath
    ) const {
        if (!std::filesystem::exists(IndexPath)) {
            return 0;
        }
        const size_t file_size = std::filesystem::file_size(IndexPath);
        if (IndexPath.extension() == kIndexFile32Ext) {
            return file_size / sizeof(uint32_t);
        } else {
            return file_size / sizeof(uint64_t);
        }
    }

    template<typename U>
    void SortIndexInternal(
        const std::filesystem::path& IndexPath,
        const std::filesystem::path& DataPath
    ) const {
        // We sort an index file by reading the whole contents into memory,
        // sorting the entries based on the corresponding products in the data file.
        // First read the index file
        std::ifstream indexFile(IndexPath, std::ios::binary);
        if (!indexFile.is_open()) {
            throw std::runtime_error("Failed to open FactorDB index file for reading: " + IndexPath.string());
        }
        const size_t indexFileSize = std::filesystem::file_size(IndexPath);
        const size_t numEntries = std::filesystem::file_size(IndexPath) / sizeof(U);
        std::vector<U> offsets(numEntries);
        indexFile.read(reinterpret_cast<char*>(offsets.data()), indexFileSize);
        indexFile.close();
        std::ifstream dataFile(DataPath, std::ios::binary);
        if (!dataFile.is_open()) {
            throw std::runtime_error("Failed to open FactorDB data file for reading: " + DataPath.string());
        }

        auto readProductAtOffset = [&](const U offset) {
            dataFile.clear();
            dataFile.seekg(static_cast<std::streamoff>(offset));
            T product = 0;
            uint8_t byte = 0;
            size_t shift = 0;
            do {
                dataFile.read(reinterpret_cast<char*>(&byte), 1);
                product |= static_cast<T>(byte & 0x7F) << shift;
                shift += 7;
            } while (byte & 0x80);
            return product;
        };

        // Comparator: sort offsets by their corresponding product in the data file.
        auto compare = [&](const U a, const U b) {
            const T productA = readProductAtOffset(a);
            const T productB = readProductAtOffset(b);
            return productA < productB;
        };
        // Sort the offsets based on the products
        std::sort(offsets.begin(), offsets.end(), compare);
        dataFile.close();
        // Now write back the sorted offsets to the index file
        std::ofstream sortedIndexFile(IndexPath, std::ios::binary | std::ios::trunc);
        if (!sortedIndexFile.is_open()) {
            throw std::runtime_error("Failed to open FactorDB index file for writing: " + IndexPath.string());
        }
        sortedIndexFile.write(reinterpret_cast<const char*>(offsets.data()), indexFileSize);
        sortedIndexFile.close();
    }

    const bool
    SortIndex(
        const std::filesystem::path& IndexPath,
        const std::filesystem::path& DataPath
    ) const {
        const size_t bits = GetIndexBitsFromPath(IndexPath);
        if (bits == 32) {
            SortIndexInternal<uint32_t>(IndexPath, DataPath);
        } else {
            SortIndexInternal<uint64_t>(IndexPath, DataPath);
        }
        return true;
    }

    static const
    std::string_view
    GetIndexExtensionForBits(
        const size_t Bits
    ) {
        if (Bits == 64) {
            return kIndexFile64Ext;
        } else {
            return kIndexFile32Ext;
        }
    }

    std::filesystem::path m_CachePath;
};

} // namespace primetools

#endif // PRIMETOOLS_FACTORDB_HPP