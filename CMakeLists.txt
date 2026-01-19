cmake_minimum_required(VERSION 3.28)

project(primetools VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(FetchContent)
FetchContent_Declare(
  googletest
  # Specify the commit you depend on and update it regularly.
  URL https://github.com/google/googletest/releases/download/v1.17.0/googletest-1.17.0.tar.gz
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

set(CMAKE_CXX_COMPILER clang++ CACHE STRING "C++ compiler" FORCE)
# set(CMAKE_CXX_COMPILER /usr/bin/clang++)
set(CMAKE_C_COMPILER clang CACHE STRING "C compiler" FORCE)
# set(CMAKE_C_COMPILER /usr/bin/clang)

set(CMAKE_CXX_FLAGS "-Wall -Werror -Wno-unused-function -Wno-deprecated-declarations -Wno-c23-extensions")
if(!APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunsafe-buffer-usage -Wunsafe-buffer-usage-in-container -mavx512f -mavx512dq")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-exceptions")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -funroll-loops -ffast-math -fomit-frame-pointer -fno-rtti")
set(CMAKE_CXX_FLAGS_DEBUG "-g -ggdb -O0")
set(CMAKE_EXE_LINKER_FLAGS "-Wc23-extensions")
if(APPLE)
    execute_process(COMMAND brew --prefix OUTPUT_VARIABLE HOMEBREW_PREFIX OUTPUT_STRIP_TRAILING_WHITESPACE)
    message("Homebrew prefix: ${HOMEBREW_PREFIX}")
    include_directories(${HOMEBREW_PREFIX}/include)
    link_directories(${HOMEBREW_PREFIX}/lib)
endif()

file(GLOB_RECURSE LIB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/lib/*.cpp)
add_library(libprimetools STATIC ${LIB_SOURCES})
target_include_directories(libprimetools PRIVATE
                                ${CMAKE_CURRENT_SOURCE_DIR}/lib/
                                ${CMAKE_CURRENT_SOURCE_DIR}/include/)

SET(PRIMETOOLS_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/primetools.cpp)
add_executable(primetools ${PRIMETOOLS_SOURCES})
target_include_directories(primetools PRIVATE
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/
                            ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_link_libraries(primetools PRIVATE libprimetools gmp gmpxx)
add_dependencies(primetools libprimetools)

SET(FACTORDB_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/factordb.cpp)
add_executable(factordb ${FACTORDB_SOURCES})
target_include_directories(factordb PRIVATE
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/
                            ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_link_libraries(factordb PRIVATE libprimetools gmp gmpxx)
add_dependencies(factordb libprimetools)

SET(ALIQUOT_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/aliquot.cpp)
add_executable(aliquot ${ALIQUOT_SOURCES})
target_include_directories(aliquot PRIVATE
                            ${CMAKE_CURRENT_SOURCE_DIR}/src/
                            ${CMAKE_CURRENT_SOURCE_DIR}/include/)
target_link_libraries(aliquot PRIVATE libprimetools gmp gmpxx)
add_dependencies(aliquot libprimetools)

# Tests
file(GLOB_RECURSE TEST_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/test/*.cpp)
add_executable(unittest EXCLUDE_FROM_ALL ${TEST_SOURCES})
target_include_directories(unittest PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/test/
    ${CMAKE_CURRENT_SOURCE_DIR}/include/
)
target_link_libraries(unittest gtest_main libprimetools gmp gmpxx)
add_test(NAME unittest COMMAND unittest)
add_dependencies(unittest libprimetools)