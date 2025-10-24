
#ifndef FIELD_H
#define FIELD_H

#include <cstddef>
#include <cstdint>

#define FIELD_WORDS 4

// Function Declarations
// Inline definitions for device code compatibility
inline __host__ __device__ void field_add(uint64_t *sum, const uint64_t *a, const uint64_t *b, size_t field_words) {
    for (size_t i = 0; i < field_words; i++)
        sum[i] = a[i] ^ b[i];
}

inline __host__ __device__ void field_addEqual(uint64_t *a, const uint64_t *b, size_t field_words) {
    for (size_t i = 0; i < field_words; i++)
        a[i] ^= b[i];
}

inline __host__ __device__ void field_sub(uint64_t *difference, const uint64_t *a, const uint64_t *b, size_t field_words) {
    field_add(difference, a, b, field_words);
}

inline __host__ __device__ void field_subEqual(uint64_t *a, const uint64_t *b, size_t field_words) {
    field_addEqual(a, b, field_words);
}

inline __host__ __device__ int device_memcmp(const uint64_t *a, const uint64_t *b, size_t num_elements) {
    for (size_t i = 0; i < num_elements; ++i) {
        if (a[i] != b[i]) {
            return (a[i] < b[i]) ? -1 : 1;
        }
    }
    return 0;
}

inline __host__ __device__ void field_mul(uint64_t *product, const uint64_t *a, const uint64_t *b, size_t field_words) {
	// This is a simplified version; for full implementation, copy from field.cu
	for (size_t i = 0; i < field_words; i++) {
		product[i] = a[i] * b[i]; // Replace with your actual field multiplication logic
	}
}

inline __host__ __device__ void field_mulEqual(uint64_t *a, const uint64_t *b, size_t field_words) {
	uint64_t temp[5];
	field_mul(temp, a, b, field_words);
	for (size_t i = 0; i < field_words; i++) {
		a[i] = temp[i];
	}
}

inline __host__ __device__ void field_one(uint64_t *one, size_t field_words) {
	for (size_t i = 0; i < field_words - 1; i++) {
		one[i] = 0;
	}
	one[field_words - 1] = 1;
}

// __host__ __device__ void field_mul(uint64_t *product, const uint64_t *a, const uint64_t *b,size_t field_words);
// __host__ __device__ void field_mulEqual(uint64_t *a, const uint64_t *b,size_t field_words);
// __host__ __device__ void field_swap_with_tmp(uint64_t *a, uint64_t *b, uint64_t *tmp, const size_t field_bytesize);
// __host__ __device__ void field_swap(uint64_t *a, uint64_t *b, const size_t field_bytesize);
// __host__ __device__ void field_one(uint64_t *one,size_t field_words);
// __host__ __device__ void field_neg(uint64_t *result, const uint64_t *a,size_t field_words);
// __host__ __device__ void all_subset_sums(uint64_t *sums, const uint64_t *elements, int /*replaced from size_t*/ num_elements, const uint64_t *shift,size_t field_words);
//__host__ __device__ uint64_t* sample(uint64_t *result, const uint8_t *byte_array, int /*replaced from size_t*/ byte_array_len, const uint64_t *basis);

// #ifdef __cplusplus
// }
// #endif


#endif // FIELD_H