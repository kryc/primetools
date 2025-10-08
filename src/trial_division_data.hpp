#include <cstdint>
#include <span>

namespace primetools {

// 30‑wheel gaps packed 16×4‑bit per 32‑bit word
static const uint32_t WHEEL30GAPSUINT32 = 0x26424246;
static const uint64_t WHEEL30GAPSDATA[1] = {
    0x0000000026424246ULL,
};
static const std::span<const uint64_t, 1> WHEEL30GAPS(WHEEL30GAPSDATA, 1);

// Wheel modulus = 210, total gaps = 48, bits_required = 4, gaps_per_word = 16, word_count = 3
alignas(uint64_t) static const uint8_t WHEEL210GAPSDATA_RAW[] = {
#embed "../rsrc/wheel210.bin"
};
static_assert(sizeof(WHEEL210GAPSDATA_RAW) == 3 * sizeof(uint64_t));
static const uint64_t* WHEEL210GAPSDATA = reinterpret_cast<const uint64_t*>(WHEEL210GAPSDATA_RAW);
static const std::span<const uint64_t, 3> WHEEL210GAPS(WHEEL210GAPSDATA, 3);

// Wheel modulus = 2310, total gaps = 480, bits_required = 4, gaps_per_word = 16, word_count = 30
alignas(uint64_t) static const uint8_t WHEEL2310GAPSDATA_RAW[] = {
#embed "../rsrc/wheel2310.bin"
};
static_assert(sizeof(WHEEL2310GAPSDATA_RAW) == 30 * sizeof(uint64_t));
static const uint64_t* WHEEL2310GAPSDATA = reinterpret_cast<const uint64_t*>(WHEEL2310GAPSDATA_RAW);
static const std::span<const uint64_t, 30> WHEEL2310GAPS(WHEEL2310GAPSDATA, 30);

// Wheel modulus = 30030, total gaps = 5760, bits_required = 4, gaps_per_word = 15, word_count = 384
alignas(uint64_t) static const uint8_t WHEEL30030GAPSDATA_RAW[] = {
#embed "../rsrc/wheel30030.bin"
};
static_assert(sizeof(WHEEL30030GAPSDATA_RAW) == 384 * sizeof(uint64_t));
static const uint64_t* WHEEL30030GAPSDATA = reinterpret_cast<const uint64_t*>(WHEEL30030GAPSDATA_RAW);
static const std::span<const uint64_t, 384> WHEEL30030GAPS(WHEEL30030GAPSDATA, 384);

// Wheel modulus = 510510, total gaps = 92160, bits_required = 4, gaps_per_word = 15, word_count = 6144
alignas(uint64_t) static const uint8_t WHEEL510510GAPSDATA_RAW[] = {
#embed "../rsrc/wheel510510.bin"
};
static_assert(sizeof(WHEEL510510GAPSDATA_RAW) == 6144 * sizeof(uint64_t));
static const uint64_t* WHEEL510510GAPSDATA = reinterpret_cast<const uint64_t*>(WHEEL510510GAPSDATA_RAW);
static const std::span<const uint64_t, 6144> WHEEL510510GAPS(WHEEL510510GAPSDATA, 6144);

// Wheel modulus = 9699690, total gaps = 1658880, bits_required = 5, gaps_per_word = 25, word_count = 66356
alignas (__uint128_t) static const int8_t WHEEL9699690GAPSDATA_RAW[] = {
#embed "../rsrc/wheel9699690.bin"
};
static_assert(sizeof(WHEEL9699690GAPSDATA_RAW) == 66356 * sizeof(__uint128_t));
static const __uint128_t* WHEEL9699690GAPSDATA = reinterpret_cast<const __uint128_t*>(WHEEL9699690GAPSDATA_RAW);
static const std::span<const __uint128_t, 66356> WHEEL9699690GAPS(reinterpret_cast<const __uint128_t*>(WHEEL9699690GAPSDATA), 66356);

// Wheel modulus = 223092870, total gaps = 36495360, bits_required = 5, gaps_per_word = 25, word_count = 1459815
alignas (__uint128_t) static const int8_t WHEEL223092870GAPSDATA_RAW[] = {
#embed "../rsrc/wheel223092870.bin"
};
static_assert(sizeof(WHEEL223092870GAPSDATA_RAW) == 1459815 * sizeof(__uint128_t));
static const __uint128_t* WHEEL223092870GAPSDATA = reinterpret_cast<const __uint128_t*>(WHEEL223092870GAPSDATA_RAW);
static const std::span<const __uint128_t, 1459815> WHEEL223092870GAPS(reinterpret_cast<const __uint128_t*>(WHEEL223092870GAPSDATA), 1459815);

} // namespace primetools