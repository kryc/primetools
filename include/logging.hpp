#ifndef LOGGING_HPP
#define LOGGING_HPP

#include <functional>
#include <iostream>
#include <string_view>

namespace primetools {

// Define a log callback type
using LogCallback = std::function<void(std::string_view)>;

static inline void LogStdOut(
    std::string_view Message
)
{
    std::cout << Message << std::endl;
}

static inline void LogQuiet(
    std::string_view Message
) { return; };

} // namespace primetools

#endif // LOGGING_HPP