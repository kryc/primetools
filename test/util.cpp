
#include <gtest/gtest.h>
#include <gmpxx.h>

#include "util.hpp"

namespace primetools {

TEST(Util, ToHex) {
    uint32_t value = 0x12345678;
    std::string hex_str = GetHex(value);
    EXPECT_EQ(hex_str, "12345678");
    value = 0xABCDEF;
    hex_str = GetHex(value);
    EXPECT_EQ(hex_str, "abcdef");
    hex_str = GetHex(value, 8);
    EXPECT_EQ(hex_str, "00abcdef");
}

} // namespace primetools