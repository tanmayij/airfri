cmake_minimum_required(VERSION 3.10)

# Benchmark executable for FRI
add_executable(benchmark_fri benchmark_fri.cpp)

target_include_directories(benchmark_fri
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/../include
    ${CMAKE_CURRENT_SOURCE_DIR}/../../additive-fft/C++
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/libff
    ${CMAKE_CURRENT_SOURCE_DIR}/../external/benchmark/include
    ${LIBSODIUM_INCLUDE_DIRS}
)

target_link_libraries(benchmark_fri
    PRIVATE
    fri_lib
    ff
    benchmark::benchmark
    ${CUDA_LIBRARIES}
    ${LIBSODIUM_LIBRARIES}
)

# Set CUDA properties
set_target_properties(benchmark_fri PROPERTIES
    CUDA_SEPARABLE_COMPILATION ON
)
