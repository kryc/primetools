#ifndef WHEEL30_HPP
#define WHEEL30_HPP

#include <cstdint>

namespace primetools {

constexpr uint32_t kWheel30 = 0x26424246;
constexpr uint32_t kWheel30BitsPerGap = 4;
constexpr uint32_t kWheel30Mask = (1 << kWheel30BitsPerGap) - 1;

} // namespace primetools

#endif // WHEEL30_HPP