# Unified Field Operations Architecture

## Overview

All field operations in this codebase now use the **proven CUDA implementation** from `field.cu`. This ensures consistency and correctness across CPU and GPU code.

## Field Representation

**Canonical Format**: `uint64_t[4]` array (little-endian)
- Word 0: bits 0-63
- Word 1: bits 64-127
- Word 2: bits 128-191
- Word 3: bits 192-255

**Field**: GF(2^256) - binary extension field

## File Structure

### Core Field Operations
- **`src/field.cu`** - Implementation (CPU + GPU compatible)
- **`include/field.cuh`** - Header with function declarations
- All functions marked `__host__ __device__` for unified compilation

### Conversion Layer
- **`include/field_conversion.hpp`** - Converts between:
  - `libff::gf256` (C++ types) ↔ `uint64_t[4]` (CUDA types)
  - Uses libff's `convert_field_element_to_bit_vector()` internally

### Usage in Different Contexts

#### 1. GPU Kernels (`commit_kernel.cu`, etc.)
```cuda
#include "../include/field.cuh"

__global__ void my_kernel(...) {
    uint64_t a[4], b[4], result[4];
    
    // use field operations directly
    field_add(result, a, b, 4);
    field_mul(result, a, b, 4);
    field_inv(result, a, 4);
}
```

#### 2. CPU Host Code with libff types
```cpp
#include <libff/algebra/fields/binary/gf256.hpp>
#include "field_conversion.hpp"
#include "field.cuh"

void host_function() {
    libff::gf256 elem1, elem2;
    
    // convert to uint64_t[4] for field operations
    uint64_t a[4], b[4], result[4];
    field_to_uint64(elem1, a);
    field_to_uint64(elem2, b);
    
    // use field operations
    field_mul(result, a, b, 4);
    
    // convert back
    libff::gf256 product = uint64_to_field(result);
}
```

#### 3. Host-Device Transfer
```cpp
// flatten vector for GPU transfer
std::vector<libff::gf256> codeword(N);
uint64_t* flat = new uint64_t[N * 4];
vector_to_uint64_array(codeword, flat);

// transfer to GPU
cudaMemcpy(d_codeword, flat, N * 4 * sizeof(uint64_t), cudaMemcpyHostToDevice);

// ... GPU kernel uses field operations ...

// transfer back
cudaMemcpy(flat, d_codeword, N * 4 * sizeof(uint64_t), cudaMemcpyDeviceToHost);
std::vector<libff::gf256> result = uint64_array_to_vector(flat, N);
```

## Available Field Operations

### Basic Arithmetic
- `field_add(result, a, b, field_words)` - Addition (XOR in GF(2^n))
- `field_sub(result, a, b, field_words)` - Subtraction (same as addition)
- `field_mul(result, a, b, field_words)` - Multiplication
- `field_neg(result, a, field_words)` - Negation (identity in GF(2^n))

### In-Place Operations
- `field_addEqual(a, b, field_words)` - a += b
- `field_subEqual(a, b, field_words)` - a -= b
- `field_mulEqual(a, b, field_words)` - a *= b

### Advanced Operations
- `field_inv(inv, base, field_words)` - Multiplicative inverse
- `field_pow(result, a, exponent, field_words)` - Exponentiation
- `field_batch_inverse_and_mul_with_precomputed(...)` - Batch inversion

### Utility Functions
- `is_zero(a, field_words)` - Check if element is zero
- `field_equal(a, b, field_words)` - Check equality
- `field_one(one, field_words)` - Get multiplicative identity
- `field_swap(a, b, field_bytesize)` - Swap two elements

## Implementation Details

### Field Multiplication (GF(2^256))
Located in `field.cu` lines 310-330 for `field_words == 4`:
- Uses `gf64_mul()` for 64-bit carryless multiplication
- Combines partial products with XOR reduction
- Proven implementation from FRI192

### Field Inversion (GF(2^256))
Located in `field.cu` lines 147-169 for `field_words == 4`:
- Computes `a^{-1} = a^{2^256 - 2}` via Fermat's Little Theorem
- Uses repeated squaring and multiplication
- Total: ~270 mul/sqr operations

## Migration Guide

### Before (inconsistent):
```cpp
// Some files used libff::gf256::operator*
FieldT result = a * b;

// Others used custom CUDA functions
field_mul_device(result, a, b);

// Conversion was unclear
```

### After (unified):
```cpp
// CPU code with libff types
#include "field_conversion.hpp"
uint64_t a_flat[4], b_flat[4], result[4];
field_to_uint64(a, a_flat);
field_to_uint64(b, b_flat);
field_mul(result, a_flat, b_flat, 4);
FieldT product = uint64_to_field(result);

// GPU code uses same field_mul
__global__ void kernel() {
    field_mul(result, a, b, 4);  // identical call!
}
```

## Current Usage

### Files Using Unified Field Operations:
✅ **`src/commit_kernel.cu`** - GPU FRI folding
- Uses: `field_add`, `field_sub`, `field_mul`, `field_inv`
- Operates on `uint64_t[4]` arrays
- Calls proven `field.cu` functions

✅ **`src/field.cu`** - Implementation
- Provides all field operations
- `__host__ __device__` compatible

✅ **`include/field_conversion.hpp`** - Conversion layer
- Bridges libff ↔ uint64_t[4]
- Used in host code before GPU transfer

### Files to Update:
⚠️ **`src/poly_eval.cpp`** - Currently uses libff operators
- Should convert to uint64_t[4] and use `field_mul` for GPU version
- CPU version can continue using libff (simpler)

⚠️ **`src/prove.cpp`** - Uses libff::gf256 operators
- Keep libff for high-level logic
- Use conversion when calling GPU kernels

## Best Practices

### Rule 1: One Source of Truth
**All field multiplication uses `field_mul()` from `field.cu`**
- Don't mix libff's `operator*` and `field_mul` in performance-critical code
- Use libff operators only for high-level host logic

### Rule 2: Boundary Conversion
**Convert at host/device boundaries, not inside kernels**
```cpp
// GOOD: Convert once before kernel
vector_to_uint64_array(codeword, flat);
my_kernel<<<...>>>(flat);

// BAD: Converting inside kernel
__global__ void my_kernel(libff::gf256* codeword) {
    // Don't do conversions here!
}
```

### Rule 3: Consistent Representation
**Always use `uint64_t[4]` in CUDA code**
- Device memory: `uint64_t*` with stride of 4
- Kernel parameters: `uint64_t*` pointers
- Local variables: `uint64_t[4]` arrays

### Rule 4: Field Words Parameter
**Always pass `FIELD_WORDS` (or 4) explicitly**
```cpp
field_mul(result, a, b, 4);  // GOOD: explicit
field_mul(result, a, b, FIELD_WORDS);  // ALSO GOOD
```

## Performance Notes

### Proven Performance from FRI192:
- Field multiplication: Optimized `gf64_mul()` with carryless multiply
- Field inversion: ~270 operations (repeated squaring)
- Memory layout: Contiguous `uint64_t[4]` for coalesced GPU access

### Optimization Opportunities:
1. **Batch operations**: Use `field_batch_inverse_and_mul_with_precomputed()`
2. **Precomputation**: Compute domain elements once, reuse
3. **Coalesced access**: Stride by `FIELD_WORDS` in device memory

## Testing Recommendations

### Validate Conversion:
```cpp
// Test: libff ↔ uint64_t[4] round-trip
libff::gf256 original = libff::gf256::random_element();
uint64_t flat[4];
field_to_uint64(original, flat);
libff::gf256 recovered = uint64_to_field(flat);
assert(original == recovered);
```

### Validate Field Operations:
```cpp
// Test: field_mul matches libff
libff::gf256 a = libff::gf256::random_element();
libff::gf256 b = libff::gf256::random_element();

// libff result
libff::gf256 libff_result = a * b;

// field.cu result
uint64_t a_flat[4], b_flat[4], result[4];
field_to_uint64(a, a_flat);
field_to_uint64(b, b_flat);
field_mul(result, a_flat, b_flat, 4);
libff::gf256 cuda_result = uint64_to_field(result);

assert(libff_result == cuda_result);
```

## Summary

✅ **Unified**: Single field implementation (`field.cu`)  
✅ **Proven**: Tested CUDA code from FRI192  
✅ **Compatible**: `__host__ __device__` works everywhere  
✅ **Clean boundaries**: libff ↔ uint64_t[4] conversion at edges  
✅ **Performance**: Optimized for GPU coalesced access  

**Bottom line**: Use `field.cu` operations everywhere. Convert libff ↔ uint64_t[4] only at boundaries.
