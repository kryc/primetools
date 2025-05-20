cmake_minimum_required(VERSION 3.28)

project(primetools VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_COMPILER clang++)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_CXX_FLAGS "-Wall -Werror -Wno-unused-function")
if(!APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunsafe-buffer-usage -Wunsafe-buffer-usage-in-container")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-exceptions")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops -ffast-math -fomit-frame-pointer -fno-rtti")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -O0")

if(APPLE)
    execute_process(COMMAND brew --prefix OUTPUT_VARIABLE HOMEBREW_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
    message("Homebrew prefix: ${HOMEBREW_PREFIX}")
    include_directories(${HOMEBREW_PREFIX}/include)
    link_directories(${HOMEBREW_PREFIX}/lib)
endif()

file(GLOB_RECURSE SOURCES "src/*.cpp")
add_executable(primetools ${SOURCES})
target_include_directories(primetools PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
target_link_libraries(primetools PRIVATE gmp gmpxx)