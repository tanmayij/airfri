#ifndef SHA3_TEST_H
#define SHA3_TEST_H

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <cuda_runtime.h>
#include "../include/hash-host.cuh" 
#include "../include/hash.cuh"    

#define BIT_SIZE 256 // Change this to test SHA3-256, 384, or 512
__device__ void SHA3(uint8_t *hm, const uint8_t *msg, size_t msg_len, size_t bitSize);
// Function to print a hash in hexadecimal format
void print_hash(const uint8_t *hash, size_t size);

// GPU kernel to compute SHA3
__global__ void sha3_kernel(const uint8_t *input, size_t input_len, uint8_t *output, size_t bitSize);

#endif // SHA3_TEST_H