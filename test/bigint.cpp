#include <array>
#include <limits>
#include <vector>

#include <gtest/gtest.h>

#include "bigint.hpp"

TEST(BigInt, Assign)
{
    BigInt<1024> bigInt;
    bigInt = 4294967295; // 2^32 - 1
    EXPECT_EQ(bigInt[0], 4294967295ull);
    EXPECT_EQ(bigInt[1], 0ull);

    mpz_class largeValue = (mpz_class(1) << 64) - 1; // 2^64 - 1
    bigInt = largeValue;
    EXPECT_EQ(bigInt[0], 0xFFFFFFFFFFFFFFFFull);
    EXPECT_EQ(bigInt[1], 0ull);

    mpz_random(largeValue.get_mpz_t(), 1024);
    bigInt = largeValue;
    for (size_t i = 0; i < 16; ++i) {
        uint64_t limb = static_cast<uint64_t>(largeValue.get_ui());
        EXPECT_EQ(bigInt[i], limb);
        largeValue >>= 64;
    }
}

TEST(BigInt, AddEquals)
{
    BigInt<1024> bigInt;
    bigInt += 1;
    EXPECT_EQ(bigInt[0], 1ull);
    bigInt = 0;
    EXPECT_EQ(bigInt[0], 0ull);
    bigInt += 0xFFFFFFFFFFFFFFFFull; // 2^64 - 1
    EXPECT_EQ(bigInt[0], 0xFFFFFFFFFFFFFFFFull);
    bigInt += 1;
    EXPECT_EQ(bigInt[0], 0ull);
    EXPECT_EQ(bigInt[1], 1ull);
    bigInt += 0xFFFFFFFFFFFFFFFFull; // 2^64 - 1
    EXPECT_EQ(bigInt[0], 0xFFFFFFFFFFFFFFFFull);
    EXPECT_EQ(bigInt[1], 1ull);
    bigInt += 0xFFFFFFFFFFFFFFFFull; // 2^64 - 1
    EXPECT_EQ(bigInt[0], 0xFFFFFFFFFFFFFFFEull);
    EXPECT_EQ(bigInt[1], 2ull);
    // Check overflow beyond 128 bits
    bigInt = 0xFFFFFFFFFFFFFFFFull;
    bigInt <<= 64;
    bigInt += 0xFFFFFFFFFFFFFFFFull;
    EXPECT_EQ(bigInt[0], 0xFFFFFFFFFFFFFFFFull);
    EXPECT_EQ(bigInt[1], 0xFFFFFFFFFFFFFFFFull);
    bigInt += 1;
    EXPECT_EQ(bigInt[0], 0u);
    EXPECT_EQ(bigInt[1], 0u);
    EXPECT_EQ(bigInt[2], 1u);
}

TEST(BigInt, SubtractEquals)
{
    BigInt<1024> bigInt;
    bigInt = 100;
    bigInt -= 50;
    EXPECT_EQ(bigInt[0], 50ull);
    bigInt -= 50;
    EXPECT_EQ(bigInt[0], 0ull);

    bigInt = 0;
    bigInt -= 1;
    EXPECT_EQ(bigInt[0], std::numeric_limits<uint64_t>::max());
    EXPECT_EQ(bigInt[1], std::numeric_limits<uint64_t>::max());

    bigInt = 0x100000000ull; // 2^32
    bigInt -= 1;
    EXPECT_EQ(bigInt[0], 0xFFFFFFFFull);
    EXPECT_EQ(bigInt[1], 0u);
}

TEST(BigInt, Divides)
{
    BigInt<1024> bigInt1, bigInt2;
    bigInt1 = 25;  // divisor
    bigInt2 = 100; // dividend
    EXPECT_TRUE(bigInt1.divides(bigInt2));
    mpz_class a, b, mult;
    // Use proper bit-size random generation (previously mpz_random used limb count, producing oversized values)
    gmp_randstate_t rs; gmp_randinit_default(rs); gmp_randseed_ui(rs, 42u);
    mpz_urandomb(a.get_mpz_t(), rs, 512); // 512-bit random
    mpz_urandomb(b.get_mpz_t(), rs, 256); // 256-bit random
    mult = a * b;
    bigInt1 = a;    // divisor
    bigInt2 = mult; // dividend
    EXPECT_TRUE(bigInt1.divides(bigInt2));
    bigInt2 += 1;
    EXPECT_FALSE(bigInt1.divides(bigInt2));
    gmp_randclear(rs);
}

TEST(BigInt, RestoringDivides)
{
    BigInt<1024> bigInt1, bigInt2;
    bigInt1 = 25;  // divisor
    bigInt2 = 100; // dividend
    EXPECT_TRUE(bigInt1.restoring_divides(bigInt2));

    gmp_randstate_t rs;
    gmp_randinit_default(rs); gmp_randseed_ui(rs, 123u);
    mpz_class a, b, mult;
    for (size_t i = 0; i < 100; ++i) {
        mpz_urandomb(a.get_mpz_t(), rs, 512);
        mpz_urandomb(b.get_mpz_t(), rs, 256);
        mult = a * b;
        bigInt1 = a;    // divisor
        bigInt2 = mult; // dividend
        EXPECT_TRUE(bigInt1.restoring_divides(bigInt2));
        bigInt2 += 1;
        EXPECT_FALSE(bigInt1.restoring_divides(bigInt2));
    }
    gmp_randclear(rs);
}