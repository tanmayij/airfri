// #include <stdio.h>
// #include <stdint.h>
// #include <string.h>
// #include "../include/hash.cuh"  // Include the GPU hash function

// // Function prototype for the SHA3 GPU implementation
// __device__ void SHA3(uint8_t *hm, const uint8_t *msg, size_t msg_len, size_t bitSize);

// __global__ void sha3_kernel(uint64_t *input1, uint64_t *input2, uint8_t *output_hash) {
//     uint8_t msg[sizeof(uint64_t) * 8];  // Buffer to store combined input

//     // Copy input1 and input2 into msg buffer
//     memcpy(msg, input1, sizeof(uint64_t) * 4);
//     memcpy(msg + sizeof(uint64_t) * 4, input2, sizeof(uint64_t) * 4);

//     // Compute SHA3-256 hash
//     SHA3(output_hash, msg, sizeof(msg), 256);
// }

// int main() {
//     // Define the two uint64_t arrays
//     uint64_t h_input1[] = {0x4905066a7907dfff, 0x9d06aab1f2646518, 0xa939797d9e8a3f82, 0x99ab3ec7c912fb69};
//     uint64_t h_input2[] = {0x4f08e48ae3b51c41, 0x6e79cc2d048334bb, 0x22575f691dcc06d6, 0x6f66b42be9d654e8};

//     uint64_t *d_input1, *d_input2;
//     uint8_t *d_hash, h_hash[32];

//     // Allocate device memory
//     cudaMalloc((void **)&d_input1, sizeof(h_input1));
//     cudaMalloc((void **)&d_input2, sizeof(h_input2));
//     cudaMalloc((void **)&d_hash, 32 * sizeof(uint8_t));

//     // Copy data from host to device
//     cudaMemcpy(d_input1, h_input1, sizeof(h_input1), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_input2, h_input2, sizeof(h_input2), cudaMemcpyHostToDevice);

//     // Launch SHA3 kernel (single thread, as SHA3 is not parallelized here)
//     sha3_kernel<<<1, 1>>>(d_input1, d_input2, d_hash);
//     cudaDeviceSynchronize();

//     // Copy the computed hash back to the host
//     cudaMemcpy(h_hash, d_hash, 32 * sizeof(uint8_t), cudaMemcpyDeviceToHost);

//     // Print the hash result
//     printf("SHA3-256 Hash: ");
//     for (size_t i = 0; i < 32; i++) {
//         printf("%02x", h_hash[i]);
//     }
//     printf("\n");

//     // Free device memory
//     cudaFree(d_input1);
//     cudaFree(d_input2);
//     cudaFree(d_hash);

//     return 0;
// }

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "../include/hash.cuh"  // Include the GPU hash function

#define HASH_WORDS 4  // SHA3-256 output size in uint64_t words

// Function prototype for the SHA3 GPU implementation
__device__ void SHA3(uint8_t *hm, const uint8_t *msg, size_t msg_len, size_t bitSize);

__global__ void sha3_kernel(uint64_t *input_value, uint8_t *output_hash, size_t input_size) {
    // Buffer to store input bytes
    uint8_t msg[16];  // Assuming 2 uint64_t inputs (16 bytes)
    memcpy(msg, input_value, input_size);

    // Compute SHA3-256 hash
    SHA3(output_hash, msg, input_size, 256);
}

// Function to check SHA3 hash on GPU
int check_sha3_hash_gpu(uint64_t *input_value, uint64_t *expected_hash, size_t input_size) {
    uint8_t *d_hash, h_hash[32];  // 32 bytes = 256-bit SHA3 output
    uint64_t *d_input;

    // Allocate device memory
    cudaMalloc((void **)&d_input, input_size);
    cudaMalloc((void **)&d_hash, 32 * sizeof(uint8_t));

    // Copy input to device
    cudaMemcpy(d_input, input_value, input_size, cudaMemcpyHostToDevice);

    // Launch SHA3 kernel (single thread, as SHA3 is not parallelized)
    sha3_kernel<<<1, 1>>>(d_input, d_hash, input_size);
    cudaDeviceSynchronize();

    // Copy hash result back to host
    cudaMemcpy(h_hash, d_hash, 32 * sizeof(uint8_t), cudaMemcpyDeviceToHost);

    // Print input bytes for debugging
    printf("GPU Input Bytes: ");
    uint8_t *byte_ptr = (uint8_t *)input_value;
    for (size_t i = 0; i < input_size; i++) {
        printf("%02x ", byte_ptr[i]);
    }
    printf("\n");

    // Print computed SHA3 hash
    printf("GPU Computed SHA3-256 Hash: ");
    for (size_t i = 0; i < 32; i++) {
        printf("%02x", h_hash[i]);
    }
    printf("\n");

    // Compare computed hash with expected hash
    int result = memcmp(expected_hash, h_hash, 32) == 0;

    // Free device memory
    cudaFree(d_input);
    cudaFree(d_hash);

    return result;
}

// Example usage
int main() {
    // Define sizes
    size_t input_size = 4 * sizeof(uint64_t);  // 4 uint64_t values (32 bytes)
    size_t hash_size = 4 * sizeof(uint64_t);   // SHA3-256 output (32 bytes)

    // Allocate memory for input values
    uint64_t *input_value = (uint64_t *)malloc(input_size);
    if (input_value == NULL) {
        fprintf(stderr, "Memory allocation failed for input_value\n");
        exit(1);
    }
    uint64_t *input_value_two = (uint64_t *)malloc(input_size);
    if (input_value_two == NULL) {
        fprintf(stderr, "Memory allocation failed for input_value\n");
        exit(1);
    }
    // Assign values to input dynamically
    //d46f02faa854a139 5eda12fdcfc6b414 4d2614da35978182 d9c7afb2bf4dc4ae
    //cc10bba2731e29a7 4191eaf10b4fbbbf dd0931c636bd24bc 35c56d2534530a85
    input_value[0] = 0x3be052336fbeb42a;
    input_value[1] = 0x955977e40235ffae;
    input_value[2] = 0x0162d69b9a4f7f8b;
    input_value[3] = 0xf322b38dcdf68013;
    input_value_two[0] = 0x3be052336fbeb42a;
    input_value_two[1] = 0x955977e40235ffae;
    input_value_two[2] = 0x0162d69b9a4f7f8b;
    input_value_two[3] = 0xf322b38dcdf68013;
    // Allocate memory for expected hash
    uint64_t *expected_hash = (uint64_t *)malloc(hash_size);
    if (expected_hash == NULL) {
        fprintf(stderr, "Memory allocation failed for expected_hash\n");
        free(input_value);
        exit(1);
    }

    // Assign expected hash values
    expected_hash[0] = 0x9b0e4a426d0ba83b;
    expected_hash[1] = 0x4ae2d3426168a81e;
    expected_hash[2] = 0x4716691a48be5757;
    expected_hash[3] = 0x1d0b857e4598b0ab;

    // Check SHA3 hash computation on GPU
    int is_match = check_sha3_hash_gpu(input_value, expected_hash, input_size);
    int is_match_too = check_sha3_hash_gpu(input_value_two, expected_hash, input_size);
    // Print the result
    if (is_match) {
        printf("GPU Hash matches the expected output!\n");
    } else {
        printf("GPU Hash does NOT match the expected output!\n");
    }

    // Free allocated memory
    free(input_value);
    free(expected_hash);
    free(input_value_two);
    return 0;
}
