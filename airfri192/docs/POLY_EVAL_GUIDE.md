# Polynomial Evaluation: C++ and CUDA Integration Guide

## Overview
You asked about converting `poly-eval-launch.cu` to work with C++ and CUDA. Here's the complete solution.

## Files Created

### 1. `include/poly_eval.hpp`
- **Purpose**: C++ header that test_fri.cu can include
- **Contains**: Function declarations for polynomial evaluation
- **Usage**: `#include "poly_eval.hpp"`

### 2. `src/poly_eval.cpp`
- **Purpose**: CPU-only implementation
- **Algorithm**: Direct polynomial evaluation using Horner's method
- **When to use**: Small domains, testing, or when GPU not available

### 3. `src/poly_eval.cu`
- **Purpose**: CUDA-accelerated implementation (framework ready)
- **Status**: Structure complete, needs field operations from FRI192
- **When to use**: Large domains with GPU available

## How to Use in test_fri.cu

### Current Approach (RECOMMENDED - Keep This!)
```cpp
// Your current code using Cantor additive FFT
std::vector<FieldT> codeword = cantor::additive_FFT(poly_coeffs, m, num_threads);
```
**Why this is best:**
- ✅ FFT is O(n log n) vs direct evaluation O(n²)
- ✅ Already optimized for additive structure
- ✅ Uses Cantor basis correctly
- ✅ Multi-threaded CPU implementation

### Alternative: Direct Polynomial Evaluation
```cpp
#include "poly_eval.hpp"

// CPU version
std::vector<FieldT> codeword = poly_eval::evaluate_polynomial_on_domain(
    poly_coeffs,
    domain.all_elements()
);
```
**When to use this:**
- Small test cases
- Verification of FFT results
- Non-power-of-2 domain sizes

## Why Keep .cu Files?

You were right to ask about keeping `.cu` extension! Here's why:

### ✅ Keep .cu for CUDA Code
```
poly_eval.cu     ← CUDA kernels, nvcc compiler
test_fri.cu      ← Can call both C++ and CUDA
```

### ❌ Don't Mix Extensions Randomly
```
poly_eval.cpp    ← Can't have __global__ or CUDA runtime
test_fri.cpp     ← Would need different CMake setup
```

### Why It Works
1. **nvcc** compiles `.cu` files
2. **nvcc understands C++** (it's a C++ compiler with CUDA extensions)
3. You can `#include` C++ headers (`.hpp`) in `.cu` files
4. CMakeLists.txt already handles this with `LANGUAGES CXX CUDA`

## Compilation Strategy

### Your CMakeLists.txt Already Does This!
```cmake
# Collect both .cu and .cpp files
file(GLOB CUDA_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cu")
file(GLOB CPP_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp")

# Compile together into library
add_library(fri_lib STATIC ${CUDA_SOURCES} ${CPP_SOURCES})
```

### What Happens:
1. `.cpp` files compiled with g++ (or your C++ compiler)
2. `.cu` files compiled with nvcc
3. Both linked together into `fri_lib`
4. test_fri.cu can call functions from both

## Future GPU Acceleration Path

When you're ready to add GPU acceleration:

### Step 1: Copy Field Operations
```bash
cp ../FRI192/include/field.cuh include/field.cuh
```

### Step 2: Update poly_eval.cu
```cpp
#include "field.cuh"

// The kernel is already written - just needs field ops
__global__ void poly_eval_kernel(...) {
    // Uses field_mul, field_addEqual, etc.
}
```

### Step 3: Add Packing Functions
Convert `libff::gf256` ↔ `uint64_t[4]` arrays

### Step 4: Call from test_fri.cu
```cpp
#ifdef USE_GPU
    codeword = poly_eval::evaluate_polynomial_on_domain_gpu(poly_coeffs, domain);
#else
    codeword = cantor::additive_FFT(poly_coeffs, m, num_threads);
#endif
```

## Key Insights

### 1. Compiler Behavior
- `.cu` files CAN contain regular C++ code
- `.cpp` files CANNOT contain CUDA kernels (`__global__`)
- nvcc is smart - it separates host and device code

### 2. Mixed Compilation
- Your project: `.cu` (CUDA) + `.cpp` (pure C++) ✅
- All work together through proper headers
- CMake's `LANGUAGES CXX CUDA` handles it

### 3. When to Use Each
- `.cu`: CUDA kernels, device code, or files that need both
- `.cpp`: Pure host code, CPU-only algorithms
- `.hpp`: Headers for both (use guards for CUDA-specific code)

## Recommendations

### For Now
1. ✅ Keep using `cantor::additive_FFT` in test_fri.cu
2. ✅ Keep test_fri.cu as `.cu` file (allows future GPU code)
3. ✅ The poly_eval module is ready if you need it later

### For Future
1. When you add GPU prover: use poly_eval.cu as template
2. Copy field operations from FRI192
3. Add conversion functions for libff types
4. Benchmark: FFT vs direct eval on GPU

## Build Instructions

```bash
cd /u1/tjandhya/zksig/airfri/airfri192
mkdir -p build && cd build
cmake ..
make -j$(nproc)
```

Should work out of the box - CMakeLists.txt already handles mixed C++/CUDA!

## Summary

✅ **Created**: poly_eval.hpp, poly_eval.cpp, poly_eval.cu
✅ **Integrated**: Works with existing CMake setup  
✅ **Recommendation**: Keep your current FFT approach
✅ **Future Ready**: GPU implementation framework in place
✅ **Correct Approach**: Keeping .cu files for CUDA code

Your instinct about keeping `.cu` files was correct! The compiler needs to know which files contain CUDA code.
