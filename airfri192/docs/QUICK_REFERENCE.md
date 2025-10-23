# Quick Reference: Unified Field Operations

## ✅ What's Complete

- **Field Operations**: `field.cu` with `field_add`, `field_sub`, `field_mul`, `field_inv`
- **Conversions**: `field_conversion.hpp` for libff ↔ uint64_t[4]
- **SHA3 Hashing**: `hash.cu` with device functions
- **Commit Kernel**: Structure complete, needs flattening code

## 📝 Cheat Sheet

### Include in CUDA Files
```cuda
#include "../include/field.cuh"          // Field operations
#include "../include/hash.cuh"           // SHA3 functions
```

### Include in C++ Files  
```cpp
#include "field_conversion.hpp"          // libff ↔ uint64_t[4]
#include <libff/algebra/fields/binary/gf256.hpp>
```

### Field Operations (CPU or GPU)
```cpp
uint64_t a[4], b[4], result[4];

field_add(result, a, b, 4);     // Addition
field_sub(result, a, b, 4);     // Subtraction
field_mul(result, a, b, 4);     // Multiplication
field_inv(result, a, 4);        // Inversion
```

### Conversion (CPU only)
```cpp
// Single element
libff::gf256 elem;
uint64_t flat[4];
field_to_uint64(elem, flat);          // libff → uint64_t[4]
libff::gf256 recovered = uint64_to_field(flat);  // uint64_t[4] → libff

// Vector
std::vector<libff::gf256> vec(N);
uint64_t* flat_array = new uint64_t[N * 4];
vector_to_uint64_array(vec, flat_array);         // batch convert
std::vector<libff::gf256> recovered = uint64_array_to_vector(flat_array, N);
```

### SHA3 (GPU)
```cuda
__device__ void my_kernel() {
    uint8_t hash_output[32];
    uint8_t* input_data;
    size_t input_len;
    
    SHA3(hash_output, input_data, input_len, 256);  // SHA3-256
}
```

## 🔧 TODO: Complete These 3 Functions

### 1. In `commit_kernel.cu` around line 220:
```cpp
// Replace placeholder with:
for(int i = 0; i < N; i++) {
    field_to_uint64(codeword[i], &flat_codeword[i * FIELD_WORDS]);
}
```

### 2. In `commit_kernel.cu` around line 240:
```cpp
// Add before flattening:
FieldT basis_inv = basis[0].inverse();
field_to_uint64(basis_inv, flat_basis_inv);
```

### 3. In `commit_kernel.cu` around line 270:
```cpp
// Replace placeholder with:
for(int i = 0; i < N/2; i++) {
    codeword_nxt[i] = uint64_to_field(&flat_codeword_nxt[i * FIELD_WORDS]);
}
```

## 📚 Documentation

- **Architecture**: `UNIFIED_FIELD_OPERATIONS.md`
- **Current Status**: `CURRENT_STATUS.md`
- **Detailed TODOs**: `COMMIT_KERNEL_TODO.md`

## ⚡ Golden Rule

**Convert at boundaries, compute with unified operations**
- Host: Use libff types → convert → pass to GPU
- GPU: Receive uint64_t[4] → compute with `field_*()` → return
- Host: Receive uint64_t[4] → convert → use as libff types
