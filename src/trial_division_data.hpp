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
static const uint64_t WHEEL210GAPSDATA[] = {
#embed "../rsrc/wheel210.bin"
};
static const std::span<const uint64_t, 3> WHEEL210GAPS(WHEEL210GAPSDATA, 3);
// Wheel modulus = 2310, total gaps = 480, bits_required = 4, gaps_per_word = 16, word_count = 30
static const uint64_t WHEEL2310GAPSDATA[] = {
#embed "../rsrc/wheel2310.bin"
};
static const std::span<const uint64_t, 30> WHEEL2310GAPS(WHEEL2310GAPSDATA, 30);
// Wheel modulus = 30030, total gaps = 5760, bits_required = 4, gaps_per_word = 15, word_count = 384
static const uint64_t WHEEL30030GAPSDATA[] = {
#embed "../rsrc/wheel30030.bin"
};
static const std::span<const uint64_t, 384> WHEEL30030GAPS(WHEEL30030GAPSDATA, 384);
// Wheel modulus = 510510, total gaps = 92160, bits_required = 4, gaps_per_word = 15, word_count = 6144
static const uint64_t WHEEL510510GAPSDATA[] = {
#embed "../rsrc/wheel510510.bin"
};
static const std::span<const uint64_t, 6144> WHEEL510510GAPS(WHEEL510510GAPSDATA, 6144);
// Wheel modulus = 9699690, total gaps = 1658880, bits_required = 5, gaps_per_word = 25, word_count = 66356
static const __int128_t WHEEL9699690GAPSDATA[] = {
#embed "../rsrc/wheel9699690.bin"
};
static const std::span<const __uint128_t, 66356> WHEEL9699690GAPS(reinterpret_cast<const __uint128_t*>(WHEEL9699690GAPSDATA), 66356);
// Wheel modulus = 223092870, total gaps = 36495360, bits_required = 5, gaps_per_word = 25, word_count = 1459815
// static const __int128_t WHEEL223092870GAPSDATA[] = {
// #embed "../rsrc/wheel223092870.bin"
// };
// static const std::span<const __uint128_t, 1459815> WHEEL223092870GAPS(reinterpret_cast<const __uint128_t*>(WHEEL223092870GAPSDATA), 1459815);

} // namespace primetools