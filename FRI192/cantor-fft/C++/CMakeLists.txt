include(CMakePrintHelpers)
cmake_minimum_required(VERSION 3.10)
project(additive-fft
    DESCRIPTION "Library for additive FFT")

# Enable C++14 or higher
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Add optimization flags based on the build type
if(CMAKE_BUILD_TYPE STREQUAL "Release")
  # Enable high-level optimizations and architecture-specific tuning
  set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -DNDEBUG")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
  # Debug flags (optional, modify as needed)
  set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -pg")
  set(CMAKE_EXE_LINKER_FLAGS "-pg")
  set(CMAKE_SHARED_LINKER_FLAGS "-pg")
endif()

cmake_print_variables(CMAKE_BUILD_TYPE)

# Find the pkg-config package
find_package(PkgConfig REQUIRED)

# Use pkg-config to find libsodium
pkg_check_modules(LIBSODIUM REQUIRED libsodium)

# Add Libff with forced Release mode
set(LIBFF_BUILD_TYPE Release CACHE STRING "" FORCE)
add_subdirectory(depends/libff)

# Add Google Benchmark with forced Release mode
set(BENCHMARK_ENABLE_GTEST_TESTS OFF CACHE BOOL "Disable GTest tests in Benchmark" FORCE)
set(BENCHMARK_ENABLE_TESTING OFF CACHE BOOL "Disable testing of the Benchmark" FORCE)
set(CMAKE_BUILD_TYPE ${LIBFF_BUILD_TYPE} CACHE STRING "" FORCE)  # Ensure dependencies use Release mode
add_subdirectory(depends/benchmark)

# Main executable
add_executable(main main.cpp)

target_include_directories(
    main

    PRIVATE
    depends/libff
    ${LIBSODIUM_INCLUDE_DIRS}
    .
)

target_link_libraries(
    main 
    PRIVATE 
    ${LIBSODIUM_LIBRARIES}
    ff
)

# Benchmark executable
add_executable(run_benchmark benchmark.cpp)

target_include_directories(
    run_benchmark

    PRIVATE
    depends/libff
    depends/benchmark/include
    ${LIBSODIUM_INCLUDE_DIRS}
    .
)

target_link_libraries(
    run_benchmark 
    PRIVATE 
    ${LIBSODIUM_LIBRARIES}
    ff
    benchmark::benchmark
)
