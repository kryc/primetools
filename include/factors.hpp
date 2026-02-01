#ifndef FACTORS_HPP
#define FACTORS_HPP

#include <algorithm>
#include <cinttypes>
#include <map>
#include <sstream>
#include <vector>

#include <gmpxx.h>

#include "maths.hpp"
#include "util.hpp"

namespace primetools {

// A placeholder for an efficient prime factor storage structure
template<typename T>
class PrimeFactors {
public:
    void AddFactor(
        const T& Factor
    ) {
        m_FactorCounts[Factor]++;
    }

    void AddFactor(
        const T& Factor,
        const size_t Count
    ) {
        m_FactorCounts[Factor] += Count;
    }

    void Update(
        const PrimeFactors<T>& Other
    ) {
        for (const auto& [prime, count] : Other.m_FactorCounts) {
            m_FactorCounts[prime] += count;
        }
    }

    static PrimeFactors<T>
    FromVector(
        const std::vector<std::pair<T, size_t>>& Factors
    ) {
        PrimeFactors<T> pf;
        for (const auto& [prime, count] : Factors) {
            pf.m_FactorCounts[prime] = count;
        }
        return pf;
    }

    static PrimeFactors<T>
    FromPair(
        const T& A,
        const T& B
    ) {
        PrimeFactors<T> pf;
        pf.AddFactor(A);
        pf.AddFactor(B);
        return pf;
    }

    static PrimeFactors<T>
    FromPair(
        const std::pair<T, T>& Pair
    ) {
        return FromPair(Pair.first, Pair.second);
    }

    bool HasFactor(
        const T& Factor
    ) const {
        return m_FactorCounts.find(Factor) != m_FactorCounts.end();
    }

    T
    LargestFactor(
        void
    ) const {
        if (m_FactorCounts.empty()) {
            return 0;
        }
        return std::max_element(
            m_FactorCounts.begin(),
            m_FactorCounts.end(),
            [](const auto& a, const auto& b) {
                return a.first < b.first;
            })->first;
    }

    size_t Count(
        void
    ) const {
        size_t total = 0;
        for (const auto& pair : m_FactorCounts) {
            total += pair.second;
        }
        return total;
    }

    size_t Size(
        void
    ) const {
        return m_FactorCounts.size();
    }

    size_t CountOf(
        const T& Factor
    ) const {
        auto it = m_FactorCounts.find(Factor);
        if (it != m_FactorCounts.end()) {
            return it->second;
        }
        return 0;
    }

    bool Empty(
        void
    ) const {
        return m_FactorCounts.empty();
    }

    void Clear(
        void
    ) {
        m_FactorCounts.clear();
    }

    template<typename T2>
    PrimeFactors<T2> Convert() const {
        PrimeFactors<T2> result;
        for (const auto& [prime, count] : m_FactorCounts) {
            result.AddFactor(primetools::ConvertType<T2>(prime), count);
        }
        return result;
    }

    std::vector<std::pair<T, size_t>> ToVector() const {
        std::vector<std::pair<T, size_t>> vec;
        for (const auto& pair : m_FactorCounts) {
            vec.emplace_back(pair.first, pair.second);
        }
        return vec;
    }

    // Serialize the factors to a byte array using VLE.
    // The format is <product><num_factors>[<factor>]+
    // We do not store the counts as these can be derived from the factor list.
    std::vector<uint8_t>
    Serialize(
        void
    ) const {
        std::vector<uint8_t> data;
        // Serialize the product
        T product = Product();
        std::vector<uint8_t> serialized = primetools::SerializeVLE(product);
        data.insert(data.end(), serialized.begin(), serialized.end());

        // Serialize the number of factors
        serialized = primetools::SerializeVLE(Count());
        data.insert(data.end(), serialized.begin(), serialized.end());

        // Serialize each factor
        for (const auto& [prime, count] : m_FactorCounts) {
            for (size_t i = 0; i < count; ++i) {
                serialized = primetools::SerializeVLE(prime);
                data.insert(data.end(), serialized.begin(), serialized.end());
            }
        }

        return data;
    }

    // Function to convert prime factors to a vector of composite factors
    // By iterating over all combinations of prime factor powers
    std::vector<T>
    GetComposite(
        const bool Sorted = false
    ) const {
        std::vector<T> composites;
        composites.push_back(1); // Start with 1 as the first composite
        for (const auto& [prime, count] : m_FactorCounts) {
            size_t current_size = composites.size();
            T prime_power = 1;
            for (size_t i = 1; i <= count; ++i) {
                prime_power *= prime;
                for (size_t j = 0; j < current_size; ++j) {
                    composites.push_back(composites[j] * prime_power);
                }
            }
        }
        if (Sorted)
        {
            std::sort(composites.begin(), composites.end());
        }
        return composites;
    }

    const T
    Product(
        void
    ) const {
        T prod = 1;
        T factor;
        for (const auto& [prime, count] : m_FactorCounts) {
            if constexpr (std::is_same_v<T, mpz_class>) {
                mpz_pow_ui(factor.get_mpz_t(), prime.get_mpz_t(), count);
            } else {
                factor = primetools::Pow(prime, count);
            }
            prod *= factor;
        }
        return prod;
    }

    uint64_t
    Product64(
        void
    ) const {
        uint64_t prod = 1;
        T factor;
        for (const auto& [prime, count] : m_FactorCounts) {
            mpz_pow_ui(factor.get_mpz_t(), prime.get_mpz_t(), count);
            prod *= factor.get_ui();
        }
        return prod;
    }

    std::string
    GetString(
        void
    ) const {
        std::stringstream oss;
        bool first = true;
        for (const auto& [prime, count] : m_FactorCounts) {
            if (!first) {
                oss << " * ";
            }
            if (count > 1)
                oss << prime << "^" << count;
            else
                oss << prime;
            first = false;
        }
        return oss.str();
    }

    const T
    MaxFactor(
        void
    ) const {
        if (m_FactorCounts.empty()) {
            return 0;
        }
        return std::max_element(
            m_FactorCounts.begin(),
            m_FactorCounts.end(),
            [](const auto& a, const auto& b) {
                return a.first < b.first;
            })->first;
    }

    size_t
    MaxPower(
        void
    ) const {
        if (m_FactorCounts.empty()) {
            return 0;
        }
        return std::max_element(
            m_FactorCounts.begin(),
            m_FactorCounts.end(),
            [](const auto& a, const auto& b) {
                return a.second < b.second;
            })->second;
    }
private:
    std::map<T, size_t> m_FactorCounts;
};

} // namespace primetools

#endif // FACTORS_HPP