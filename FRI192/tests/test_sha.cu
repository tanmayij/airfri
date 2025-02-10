#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "test_sha.cuh"

#define BIT_SIZE 256 // Change this to test SHA3-256, 384, or 512

// Function to print a hash in hexadecimal format
void print_hash(const uint8_t *hash, size_t size) {
    for (size_t i = 0; i < size; i++) {
        printf("%02x", hash[i]);
    }
    printf("\n");
}

// GPU kernel to compute SHA3
__global__ void sha3_kernel(const uint8_t *input, size_t input_len, uint8_t *output, size_t bitSize) {
    SHA3(output, input, input_len, bitSize);
}

int main() {
    // Input: uint64_t array of 4 elements
    uint64_t input[4] = {
        0x1122334455667788ULL,
        0x99aabbccddeeff00ULL,
        0x123456789abcdef0ULL,
        0xfedcba9876543210ULL
    };
    uint8_t *byte_input = (uint8_t *)input;
    size_t input_len = sizeof(input);

    // Host output buffer for SHA3_host
    uint8_t host_output[BIT_SIZE / 8];

    // Compute SHA3 on the host
    SHA3_host(host_output, byte_input, input_len, BIT_SIZE);

    // Allocate memory on the GPU for input and output
    uint8_t *d_input, *d_output;
    size_t output_len = BIT_SIZE / 8;

    cudaMalloc(&d_input, input_len);
    cudaMalloc(&d_output, output_len);

    // Copy input data to the GPU
    cudaMemcpy(d_input, byte_input, input_len, cudaMemcpyHostToDevice);

    // Launch kernel to compute SHA3 on the GPU
    sha3_kernel<<<1, 1>>>(d_input, input_len, d_output, BIT_SIZE);

    // Copy the result back to the host
    uint8_t gpu_output[BIT_SIZE / 8];
    cudaMemcpy(gpu_output, d_output, output_len, cudaMemcpyDeviceToHost);

    // Free GPU memory
    cudaFree(d_input);
    cudaFree(d_output);

    // Print and compare results
    printf("SHA3-%d Host Output: ", BIT_SIZE);
    print_hash(host_output, output_len);

    printf("SHA3-%d GPU Output:  ", BIT_SIZE);
    print_hash(gpu_output, output_len);

    // Verify if outputs match
    if (memcmp(host_output, gpu_output, output_len) == 0) {
        printf("Success: Host and GPU outputs match!\n");
    } else {
        printf("Error: Host and GPU outputs do not match!\n");
    }

    return 0;
}