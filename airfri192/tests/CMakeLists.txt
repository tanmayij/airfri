cmake_minimum_required(VERSION 3.14)

# Find required packages
find_package(Threads REQUIRED)

# Include directories
include_directories(
    ${CMAKE_SOURCE_DIR}
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/external/libff
    ${CMAKE_SOURCE_DIR}/external/libiop
    ${CMAKE_SOURCE_DIR}/../additive-fft/C++
)

# Create executable
add_executable(test_airfri test_fri.cpp)

# Link libraries
target_link_libraries(test_airfri
    Threads::Threads
    ff
    iop
)
