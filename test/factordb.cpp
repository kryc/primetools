#include <filesystem>

#include <gtest/gtest.h>
#include <gmpxx.h>

#include "factordb.hpp"

namespace primetools {

TEST(FactorDB, GetPaths) {
    FactorDB<mpz_class> db("/tmp/factordb_cache");
    db.Clear();
    mpz_class value(0x123456789abcdef);
    auto index_path = db.GetFactorIndex(value);
    EXPECT_EQ(index_path, std::filesystem::path("/tmp/factordb_cache/factors_bcdef.i32"));
    auto data_path = db.GetFactorData(value);
    EXPECT_EQ(data_path, std::filesystem::path("/tmp/factordb_cache/factors_bcdef.dat"));
    value = 0x01;
    index_path = db.GetFactorIndex(value);
    EXPECT_EQ(index_path, std::filesystem::path("/tmp/factordb_cache/factors_00001.i32"));
    data_path = db.GetFactorData(value);
    EXPECT_EQ(data_path, std::filesystem::path("/tmp/factordb_cache/factors_00001.dat"));
}

TEST(FactorDB, WriteAndReadFactors) {
    FactorDB<mpz_class> db("/tmp/factordb_cache");
    db.Clear();
    PrimeFactors<mpz_class> factors;
    factors.AddFactor(2_mpz);
    factors.AddFactor(2_mpz);
    factors.AddFactor(3_mpz);
    factors.AddFactor(5_mpz);
    factors.AddFactor(7_mpz);
    db.AddFactors(factors);
    // Product = 2 * 2 * 3 * 5 * 7 = 420

    // Make sure the files exist
    auto index_path = db.GetFactorIndex(factors.Product());
    auto data_path = db.GetFactorData(factors.Product());
    EXPECT_TRUE(std::filesystem::exists(index_path));
    EXPECT_TRUE(std::filesystem::exists(data_path));

    mpz_class product = factors.Product();
    auto read_factors_opt = db.GetFactors(product);
    ASSERT_TRUE(read_factors_opt.has_value());
    PrimeFactors<mpz_class> read_factors = read_factors_opt.value();
    EXPECT_EQ(read_factors.CountOf(2_mpz), 2);
    EXPECT_EQ(read_factors.CountOf(3_mpz), 1);
    EXPECT_EQ(read_factors.CountOf(5_mpz), 1);
    EXPECT_EQ(read_factors.CountOf(7_mpz), 1);
    EXPECT_EQ(read_factors.Count(), 5);
    EXPECT_EQ(read_factors.Size(), 4);
    EXPECT_EQ(read_factors.Product(), 420);
}

} // namespace primetools