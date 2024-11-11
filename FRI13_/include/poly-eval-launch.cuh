#ifndef POLY_EVAL_CUH
#define POLY_EVAL_CUH

#include <cuda_runtime.h>
#include <stdint.h>
#include <stdio.h>

// Kernel declaration for performing poly_eval in parallel on the GPU
__global__ void poly_eval_kernel(uint64_t *codeword, uint64_t *poly, int poly_len, uint64_t *domain_elements, size_t field_words, size_t initial_domain_length);

// CUDA wrapper function to handle parallel poly_eval on the GPU
void parallel_poly_eval(uint64_t **codeword, uint64_t *poly_coeffs, int poly_len, uint64_t **domain_elements, size_t field_words, size_t initial_domain_length);

#endif // POLY_EVAL_CUH