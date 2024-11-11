#ifndef FIELD_H
#define FIELD_H

#include <stddef.h>
#include <stdint.h>

#define FIELD_WORDS 5

// Function Declarations

#ifdef __cplusplus
extern "C" {
#endif
// int field_words = FIELD_WORDS;

// Unified CPU/GPU functions
//__host__ __device__ int device_memcmp(const uint64_t *a, const uint64_t *b, int /*replaced from size_t*/ num_elements);
__host__ __device__ int is_zero(const uint64_t *a,size_t field_words);
__host__ __device__ int field_equal(const uint64_t *a, const uint64_t *b,size_t field_words);
__host__ __device__ void field_pow(uint64_t *result, const uint64_t *a, const size_t exponent,size_t field_words);
__host__ __device__ void field_inv(uint64_t *inv, const uint64_t *base,size_t field_words);
__host__ __device__ void field_batch_inverse_and_mul_with_precomputed(uint64_t *result, const uint64_t *elements, const size_t elements_len, const uint64_t *mul,size_t field_words, const uint64_t *precomputed_inverses, int /*replaced from size_t*/ num_precomputed);
__host__ __device__ void field_add(uint64_t *sum, const uint64_t *a, const uint64_t *b,size_t field_words);
__host__ __device__ void field_addEqual(uint64_t *a, const uint64_t *b, size_t field_words);
__host__ __device__ void field_sub(uint64_t *difference, const uint64_t *a, const uint64_t *b,size_t field_words);
__host__ __device__ void field_subEqual(uint64_t *a, const uint64_t *b,size_t field_words);
__host__ __device__ void field_mul(uint64_t *product, const uint64_t *a, const uint64_t *b,size_t field_words);
__host__ __device__ void field_mulEqual(uint64_t *a, const uint64_t *b,size_t field_words);
__host__ __device__ void field_swap_with_tmp(uint64_t *a, uint64_t *b, uint64_t *tmp, const size_t field_bytesize);
__host__ __device__ void field_swap(uint64_t *a, uint64_t *b, const size_t field_bytesize);
__host__ __device__ void field_one(uint64_t *one,size_t field_words);
__host__ __device__ void field_neg(uint64_t *result, const uint64_t *a,size_t field_words);
__host__ __device__ void all_subset_sums(uint64_t *sums, const uint64_t *elements, int /*replaced from size_t*/ num_elements, const uint64_t *shift,size_t field_words);
//__host__ __device__ uint64_t* sample(uint64_t *result, const uint8_t *byte_array, int /*replaced from size_t*/ byte_array_len, const uint64_t *basis);

#ifdef __cplusplus
}
#endif

#endif // FIELD_H