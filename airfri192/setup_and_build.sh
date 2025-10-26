#!/bin/bash

# Script to set up airfri192 project dependencies and build

echo "Setting up airfri192 project..."

# Navigate to airfri192 directory
cd "$(dirname "$0")"

# Check if libff exists, if not, add it as submodule
if [ ! -d "external/libff/.git" ]; then
    echo "Adding libff as submodule..."
    git submodule add https://github.com/scipr-lab/libff.git external/libff
    git submodule update --init --recursive
else
    echo "libff submodule already exists"
fi

# Check if Google Benchmark exists
if [ ! -d "external/benchmark/.git" ]; then
    echo "Adding Google Benchmark as submodule..."
    git submodule add https://github.com/google/benchmark.git external/benchmark
    git submodule update --init --recursive
else
    echo "Google Benchmark submodule already exists"
fi

echo "Submodules setup complete!"

# Create build directory
if [ ! -d "build" ]; then
    mkdir build
    echo "Created build directory"
fi

# Build the project
echo "Building project..."
cd build
cmake ..
make test_fri

echo "Build complete!"
echo "Run the test with: ./build/tests/test_fri"
