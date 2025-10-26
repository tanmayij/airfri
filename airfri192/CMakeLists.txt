cmake_minimum_required(VERSION 3.14)
project(airfri192
    DESCRIPTION "FRI implementation for 192-bit security with CUDA"
    LANGUAGES CXX CUDA)

# Enable C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 14)
set(CMAKE_CUDA_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set build type
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Add optimization flags
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
  set(CMAKE_CUDA_FLAGS_RELEASE "-O3 -DNDEBUG")
endif()

# Find CUDA
find_package(CUDA REQUIRED)

# Find libsodium
find_package(PkgConfig REQUIRED)
pkg_check_modules(LIBSODIUM REQUIRED libsodium)

# Find GMP
find_library(GMP_LIBRARY NAMES gmp REQUIRED)
find_library(GMPXX_LIBRARY NAMES gmpxx REQUIRED)

# Disable all tests and benchmarks for dependencies
set(WITH_PROCPS OFF CACHE BOOL "" FORCE)
set(IS_LIBFF_PARENT ON CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "" FORCE)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "" FORCE)

# Add libff first (this is the main one we want)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/libff/CMakeLists.txt")
    add_subdirectory(external/libff EXCLUDE_FROM_ALL)
else()
    message(FATAL_ERROR "libff not found")
endif()

# Add libfqfft (dependency of libiop)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/external/libiop/depends/libfqfft/CMakeLists.txt")
    set(IS_LIBFQFFT_PARENT ON CACHE BOOL "" FORCE)
    add_subdirectory(external/libiop/depends/libfqfft EXCLUDE_FROM_ALL)
else()
    message(FATAL_ERROR "libfqfft not found")
endif()

# Include directories
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/external/libff
    ${CMAKE_CURRENT_SOURCE_DIR}/external/libiop
    ${CMAKE_CURRENT_SOURCE_DIR}/external/libiop/depends/libfqfft
    ${CMAKE_CURRENT_SOURCE_DIR}/../additive-fft/C++
    ${CUDA_INCLUDE_DIRS}
    ${LIBSODIUM_INCLUDE_DIRS}
)

# Explicitly list only the libiop source files we need
set(LIBIOP_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/external/libiop/libiop/algebra/utils.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/external/libiop/libiop/bcs/hashing/blake2b.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/external/libiop/libiop/common/common.cpp
)

# Create libiop library manually
add_library(iop STATIC ${LIBIOP_SOURCES})
target_include_directories(iop PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/external/libiop
    ${CMAKE_CURRENT_SOURCE_DIR}/external/libff
    ${CMAKE_CURRENT_SOURCE_DIR}/external/libiop/depends/libfqfft
)
target_link_libraries(iop PUBLIC ff ${LIBSODIUM_LIBRARIES})

# Collect source files for FRI
file(GLOB CUDA_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu")
file(GLOB CPP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

# Create FRI library if sources exist
if(CUDA_SOURCES OR CPP_SOURCES)
    add_library(fri_lib STATIC ${CUDA_SOURCES} ${CPP_SOURCES})
    target_include_directories(fri_lib PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/../additive-fft/C++
        ${LIBSODIUM_INCLUDE_DIRS}
    )
    target_link_libraries(fri_lib PUBLIC
        ${CUDA_LIBRARIES}
        ${LIBSODIUM_LIBRARIES}
        ff
    )
endif()

# Add tests subdirectory  
add_subdirectory(tests)
