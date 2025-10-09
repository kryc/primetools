#include <cstdint>
#include <span>

namespace primetools {

typedef enum _PackingType {
    Unpacked = 0,
    FastPack = 1,
    DensePack = 2
} PackingType;

// 30‑wheel gaps packed 16×4‑bit per 32‑bit word
static const uint32_t WHEEL30GAPSUINT32 = 0x26424246;
static const uint64_t WHEEL30GAPSDATA[1] = {
    0x0000000026424246ULL,
};
static constexpr size_t WHEEL30GAP_COUNT = 1;
static const std::span<const uint64_t, 1> WHEEL30GAPS(WHEEL30GAPSDATA, 1);

// Wheel modulus = 210, total gaps = 48, bits_required = 4, gaps_per_word = 16, word_count = 3
alignas(uint64_t) static const int8_t WHEEL210GAPSDATA_RAW[] = {
#embed "../rsrc/wheel210.bin"
};
static constexpr size_t WHEEL210GAP_COUNT = sizeof(WHEEL210GAPSDATA_RAW) / sizeof(uint64_t);
static const uint64_t* WHEEL210GAPSDATA = reinterpret_cast<const uint64_t*>(WHEEL210GAPSDATA_RAW);
static const std::span<const uint64_t, 3> WHEEL210GAPS(WHEEL210GAPSDATA, 3);

// Wheel modulus = 2310, total gaps = 480, bits_required = 4, gaps_per_word = 16, word_count = 30
alignas(uint64_t) static const int8_t WHEEL2310GAPSDATA_RAW[] = {
#embed "../rsrc/wheel2310.bin"
};
static constexpr size_t WHEEL2310GAP_COUNT = sizeof(WHEEL2310GAPSDATA_RAW) / sizeof(uint64_t);
static const uint64_t* WHEEL2310GAPSDATA = reinterpret_cast<const uint64_t*>(WHEEL2310GAPSDATA_RAW);
static const std::span<const uint64_t, WHEEL2310GAP_COUNT> WHEEL2310GAPS(WHEEL2310GAPSDATA, WHEEL2310GAP_COUNT);

// Wheel modulus = 30030, total gaps = 5760, bits_required = 4, gaps_per_word = 15, word_count = 360
alignas(uint64_t) static const int8_t WHEEL30030GAPSDATA_RAW[] = {
#embed "../rsrc/wheel30030.bin"
};
static constexpr size_t WHEEL30030GAP_COUNT = sizeof(WHEEL30030GAPSDATA_RAW) / sizeof(uint64_t);
static const uint64_t* WHEEL30030GAPSDATA = reinterpret_cast<const uint64_t*>(WHEEL30030GAPSDATA_RAW);
static const std::span<const uint64_t, WHEEL30030GAP_COUNT> WHEEL30030GAPS(WHEEL30030GAPSDATA, WHEEL30030GAP_COUNT);

// Wheel modulus = 510510, total gaps = 92160, bits_required = 4, gaps_per_word = 16, word_count = 5760
alignas(uint64_t) static const int8_t WHEEL510510GAPSDATA_RAW[] = {
#embed "../rsrc/wheel510510.bin"
};
static constexpr size_t WHEEL510510GAP_COUNT = sizeof(WHEEL510510GAPSDATA_RAW) / sizeof(uint64_t);
static const uint64_t* WHEEL510510GAPSDATA = reinterpret_cast<const uint64_t*>(WHEEL510510GAPSDATA_RAW);
static const std::span<const uint64_t, WHEEL510510GAP_COUNT> WHEEL510510GAPS(WHEEL510510GAPSDATA, WHEEL510510GAP_COUNT);

// Wheel modulus = 9699690, total gaps = 1658880, bits_required = 5, gaps_per_word = 25, word_count = 66356
alignas (__uint128_t) static const int8_t WHEEL9699690GAPSDATA_RAW[] = {
#embed "../rsrc/wheel9699690.bin"
};
static constexpr size_t WHEEL9699690GAP_COUNT = sizeof(WHEEL9699690GAPSDATA_RAW) / sizeof(__uint128_t);
static const __uint128_t* WHEEL9699690GAPSDATA = reinterpret_cast<const __uint128_t*>(WHEEL9699690GAPSDATA_RAW);
static const std::span<const __uint128_t, WHEEL9699690GAP_COUNT> WHEEL9699690GAPS(reinterpret_cast<const __uint128_t*>(WHEEL9699690GAPSDATA), WHEEL9699690GAP_COUNT);

// Wheel modulus = 223092870, total gaps = 36495360, bits_required = 5, gaps_per_word = 25, word_count = 1459815
alignas (__uint128_t) static const int8_t WHEEL223092870GAPSDATA_RAW[] = {
#embed "../rsrc/wheel223092870.bin"
};
static constexpr size_t WHEEL223092870GAP_COUNT = sizeof(WHEEL223092870GAPSDATA_RAW) / sizeof(__uint128_t);
static const __uint128_t* WHEEL223092870GAPSDATA = reinterpret_cast<const __uint128_t*>(WHEEL223092870GAPSDATA_RAW);
static const std::span<const __uint128_t, WHEEL223092870GAP_COUNT> WHEEL223092870GAPS(reinterpret_cast<const __uint128_t*>(WHEEL223092870GAPSDATA), WHEEL223092870GAP_COUNT);

} // namespace primetools