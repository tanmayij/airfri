#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "../include/hash.cuh"  // Include the GPU hash function

#define HASH_WORDS 4  // SHA3-256 output size in uint64_t words

// Function prototype for the SHA3 GPU implementation
__device__ void SHA3(uint8_t *hm, const uint8_t *msg, size_t msg_len, size_t bitSize);

__global__ void sha3_kernel(uint64_t *input_value, uint64_t *output_hash, size_t input_size) {
    uint8_t msg[16];  // Assuming 2 uint64_t inputs (16 bytes)
    memcpy(msg, input_value, input_size);

    // Compute SHA3-256 hash
    SHA3((uint8_t *)output_hash, msg, input_size, 256);
}

// Function to check SHA3 hash on GPU
int check_sha3_hash_gpu(uint64_t *input_value, uint64_t *expected_hash, size_t input_size) {
    uint64_t *d_hash, h_hash[HASH_WORDS];  // SHA3-256 output (4 uint64_t values)
    uint64_t *d_input;

    // Allocate device memory
    cudaMalloc((void **)&d_input, input_size);
    cudaMalloc((void **)&d_hash, HASH_WORDS * sizeof(uint64_t));

    // Copy input to device
    cudaMemcpy(d_input, input_value, input_size, cudaMemcpyHostToDevice);

    // Launch SHA3 kernel (single thread, as SHA3 is not parallelized)
    sha3_kernel<<<1, 1>>>(d_input, d_hash, input_size);
    cudaDeviceSynchronize();

    // Copy hash result back to host
    cudaMemcpy(h_hash, d_hash, HASH_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);

    // Print input values for debugging
    printf("GPU Input Values: ");
    for (size_t i = 0; i < input_size / sizeof(uint64_t); i++) {
        printf("%016lx ", input_value[i]);
    }
    printf("\n");

    // Print computed SHA3 hash as uint64_t
    printf("GPU Computed SHA3-256 Hash (as uint64_t): ");
    for (size_t i = 0; i < HASH_WORDS; i++) {
        printf("%016lx ", h_hash[i]);
    }
    printf("\n");

    // Compare computed hash with expected hash
    int result = memcmp(expected_hash, h_hash, HASH_WORDS * sizeof(uint64_t)) == 0;

    // Free device memory
    cudaFree(d_input);
    cudaFree(d_hash);

    return result;
}

// Example usage
int main() {
    // Define sizes
    size_t input_size = 16 * sizeof(uint64_t);  // 16 uint64_t values (128 bytes)
    size_t hash_size = 4 * sizeof(uint64_t);    // SHA3-256 output (32 bytes)

    // Allocate memory for input values
    uint64_t *input_value = (uint64_t *)malloc(input_size);
    if (input_value == NULL) {
        fprintf(stderr, "Memory allocation failed for input_value\n");
        exit(1);
    }

    // Assign values to input dynamically

    input_value[0] = 0xea09e462003c9a76;
    input_value[1] = 0x7f6202591b6d54ab;
    input_value[2] = 0xb5f1ebb26fe44a2f;
    input_value[3] = 0x5701e630a806778c;
    input_value[4] = 0x1ecd7a4ab500524e;
    input_value[5] = 0x890005156c5daab9;
    input_value[6] = 0x484ee268faa0f3af;
    input_value[7] = 0x732f3b116510ec16;
    input_value[8] = 0xb5482ae7600392bf;
    input_value[9] = 0x255dbd5446ab15c8;
    input_value[10] = 0x3f1476e9a98a36f1;
    input_value[11] = 0x364600387d9e9355;
    input_value[12] = 0x40f4af00a3929a87;
    input_value[13] = 0xe3044a26718c63a2;
    input_value[14] = 0x7d3e2f4460b1628e;
    input_value[15] = 0x8d74c72527cc8b27;

    // Allocate memory for expected hash
    uint64_t *expected_hash = (uint64_t *)malloc(hash_size);
    if (expected_hash == NULL) {
        fprintf(stderr, "Memory allocation failed for expected_hash\n");
        free(input_value);
        exit(1);
    }

    // Assign expected hash values
    expected_hash[0] = 0x807ef22d564090aa;
    expected_hash[1] = 0x192c3385ff3473df;
    expected_hash[2] = 0xc54461a269b5e70b;
    expected_hash[3] = 0x64cb44f4f0e32af0;

    // Check SHA3 hash computation on GPU
    int is_match = check_sha3_hash_gpu(input_value, expected_hash, input_size);

    // Print the result
    if (is_match) {
        printf("GPU Hash matches the expected output!\n");
    } else {
        printf("GPU Hash does NOT match the expected output!\n");
    }

    // Free allocated memory
    free(input_value);
    free(expected_hash);

    return 0;
}