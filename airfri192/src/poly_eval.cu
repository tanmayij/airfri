#include <cuda_runtime.h>
#include <stdint.h>
#include <stdio.h>
#include <cstring>
#include "poly_eval.hpp"

#define FIELD_WORDS 4

// Forward declaration of field operations (you'll need to include or implement these)
__device__ void field_mul(uint64_t *result, const uint64_t *a, const uint64_t *b, size_t field_words);
__device__ void field_mulEqual(uint64_t *a, const uint64_t *b, size_t field_words);
__device__ void field_addEqual(uint64_t *a, const uint64_t *b, size_t field_words);

// Kernel to perform poly_eval in parallel on the GPU
__global__ void poly_eval_kernel(
    uint64_t *codeword,        // Flattened 1D array
    uint64_t *poly, 
    int poly_len, 
    uint64_t *domain_elements, // Flattened 1D array
    size_t field_words, 
    size_t initial_domain_length
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < initial_domain_length) {
        // Compute the starting index for the i-th element in the flattened codeword array
        uint64_t *result = &codeword[i * field_words];

        // Initialize the result to zero
        for (size_t k = 0; k < field_words; k++) {
            result[k] = 0;
        }

        uint64_t term[FIELD_WORDS];  // Temporary array for term calculation

        // Perform polynomial evaluation
        for (size_t j = 0; j < poly_len; j++) {
            if (j == 0) {
                // Copy poly[j] to term
                for (size_t k = 0; k < field_words; k++) {
                    term[k] = poly[j * field_words + k];
                }
            } else {
                // Compute the starting index for the i-th domain element
                uint64_t *domain_elem = &domain_elements[i * field_words];

                // Multiply term by poly[j] and domain_elements[i]
                field_mul(term, &poly[j * field_words], domain_elem, field_words);

                // Repeatedly multiply by domain_elements[i] if necessary
                for (size_t k = 1; k < j; k++) {
                    field_mulEqual(term, domain_elem, field_words);
                }
            }

            // Accumulate the result
            field_addEqual(result, term, field_words);
        }
    }
}

// CUDA wrapper function - C linkage for calling from C++
extern "C" void parallel_poly_eval_cuda(
    uint64_t **codeword, 
    uint64_t *poly_coeffs, 
    int poly_len, 
    uint64_t **domain_elements, 
    size_t field_words, 
    size_t initial_domain_length
) {
    // GPU memory allocations
    uint64_t *d_codeword, *d_poly_coeffs, *d_domain_elements;

    // Flatten the codeword and domain_elements on the host
    uint64_t *flattened_codeword = (uint64_t *)malloc(initial_domain_length * field_words * sizeof(uint64_t));
    uint64_t *flattened_domain_elements = (uint64_t *)malloc(initial_domain_length * field_words * sizeof(uint64_t));

    // Flatten codeword
    for (size_t i = 0; i < initial_domain_length; ++i) {
        for (size_t j = 0; j < field_words; ++j) {
            flattened_codeword[i * field_words + j] = codeword[i][j];
        }
    }

    // Flatten domain_elements
    for (size_t i = 0; i < initial_domain_length; ++i) {
        for (size_t j = 0; j < field_words; ++j) {
            flattened_domain_elements[i * field_words + j] = domain_elements[i][j];
        }
    }

    // Allocate memory on the GPU for codeword, poly_coeffs, and domain_elements
    cudaMalloc((void **)&d_codeword, initial_domain_length * field_words * sizeof(uint64_t));
    cudaMalloc((void **)&d_poly_coeffs, poly_len * field_words * sizeof(uint64_t));
    cudaMalloc((void **)&d_domain_elements, initial_domain_length * field_words * sizeof(uint64_t));

    // Copy data from host to device
    cudaMemcpy(d_codeword, flattened_codeword, initial_domain_length * field_words * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_poly_coeffs, poly_coeffs, poly_len * field_words * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_domain_elements, flattened_domain_elements, initial_domain_length * field_words * sizeof(uint64_t), cudaMemcpyHostToDevice);

    // Define the number of threads and blocks
    int threads_per_block = 256;
    int num_blocks = ((initial_domain_length + threads_per_block - 1) / threads_per_block) * 2; //512 because 512 * 256 = 131072

    // Launch the kernel
    poly_eval_kernel<<<num_blocks, threads_per_block>>>(d_codeword, d_poly_coeffs, poly_len, d_domain_elements, field_words, initial_domain_length);

    // Wait for the kernel to finish
    cudaDeviceSynchronize();

    // Copy the results back to the host (unflatten)
    cudaMemcpy(flattened_codeword, d_codeword, initial_domain_length * field_words * sizeof(uint64_t), cudaMemcpyDeviceToHost);

    // Unflatten the codeword back into the original 2D array
    for (size_t i = 0; i < initial_domain_length; ++i) {
        for (size_t j = 0; j < field_words; ++j) {
            codeword[i][j] = flattened_codeword[i * field_words + j];
        }
    }

    // Free GPU memory
    cudaFree(d_codeword);
    cudaFree(d_poly_coeffs);
    cudaFree(d_domain_elements);

    // Free host flattened arrays
    free(flattened_codeword);
    free(flattened_domain_elements);
}

// C++ wrapper that converts libff::gf256 vectors to raw arrays
namespace poly_eval {

std::vector<libff::gf256> evaluate_polynomial_on_domain_gpu(
    const std::vector<libff::gf256>& poly_coeffs,
    const std::vector<libff::gf256>& domain_elements
) {
    // This is a placeholder for GPU-accelerated version
    // You would need to:
    // 1. Pack libff::gf256 vectors to uint64_t arrays
    // 2. Call parallel_poly_eval_cuda
    // 3. Unpack results back to libff::gf256
    
    printf("GPU polynomial evaluation not yet implemented with libff types\n");
    printf("Using CPU fallback\n");
    
    // For now, fall back to CPU version
    return evaluate_polynomial_on_domain(poly_coeffs, domain_elements);
}

} // namespace poly_eval
