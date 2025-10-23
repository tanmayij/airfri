# Example usage in test_fri.cu showing how to integrate poly_eval

## Option 1: Use CPU version (already created)
Just include the header and use the function:

```cpp
#include "poly_eval.hpp"

// In test_fri.cu:
std::vector<FieldT> codeword = poly_eval::evaluate_polynomial_on_domain(
    poly_coeffs, 
    domain.all_elements()
);
```

## Option 2: Use your existing Cantor FFT (recommended for now)
Since you already have the additive FFT working with Cantor basis:

```cpp
std::vector<FieldT> codeword = cantor::additive_FFT(poly_coeffs, m, num_threads);
```

## Option 3: Later, use GPU-accelerated poly eval
When you're ready to add GPU acceleration:

1. Copy field operations from FRI192/include/field.cuh to airfri192/include/
2. Update poly_eval.cu to include field operations
3. Call the GPU version:

```cpp
#ifdef USE_GPU_POLY_EVAL
    std::vector<FieldT> codeword = poly_eval::evaluate_polynomial_on_domain_gpu(
        poly_coeffs,
        domain.all_elements()
    );
#else
    std::vector<FieldT> codeword = cantor::additive_FFT(poly_coeffs, m, num_threads);
#endif
```

## Files created:
- include/poly_eval.hpp - C++ interface header
- src/poly_eval.cpp - CPU implementation
- src/poly_eval.cu - CUDA implementation (for future GPU acceleration)

## Next steps:
1. For now, keep using cantor::additive_FFT - it's working and efficient
2. When you want GPU acceleration, we'll need to:
   - Copy field operations from FRI192 to airfri192
   - Add conversion functions for libff::gf256 ↔ uint64_t arrays
   - Update CMakeLists.txt to compile the .cu file with nvcc
