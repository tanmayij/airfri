#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <cmath>
#include "../include/hash.cuh"
#include "../include/field_conversion.hpp"
#include "../include/field.cuh"
#include <libff/algebra/fields/binary/gf256.hpp>

#define FIELD_WORDS 4
#define HASH_WORDS 4
#define CONCAT_WORDS 8  //field + hash

using FieldT = libff::gf256;

// use proven field operations from field.cu
// these are already __host__ __device__ compatible

//compute i-th element in span of cantor basis
__device__ void i_th_element_in_span(uint64_t* result, const uint64_t* basis, int basis_len, int i) {
    //initialize to zero
    for(int w = 0; w < FIELD_WORDS; w++) {
        result[w] = 0;
    }
    
    //add basis[j] for each bit j set in i
    for(int bit = 0; bit < basis_len; bit++) {
        if(i & (1 << bit)) {
            //basis is flattened, so basis[bit] is at basis[bit * FIELD_WORDS]
            for(int w = 0; w < FIELD_WORDS; w++) {
                result[w] ^= basis[bit * FIELD_WORDS + w];
            }
        }
    }
}

//helper for debugging
__device__ void print_field_device(const char* label, const uint64_t* field) {
    printf("%s: ", label);
    for(int i = 0; i < FIELD_WORDS; i++) {
        printf("%016lx ", field[i]);
    }
    printf("\n");
}

//main commit kernel - computes folded codeword
//cantor basis simplification: next[j] = current[2j] + alpha * (current[2j+1] - current[2j])
//the denominator (alpha - domain[2j]) cancels out in the cantor basis structure
__global__ void commit_kernel(
    uint64_t* device_codeword,          //input: current codeword (flattened)
    uint64_t* device_codeword_nxt,      //output: next codeword (folded)
    uint64_t* device_alpha,             //challenge field element
    uint64_t* device_offset,            //affine shift (unused in simplified formula)
    uint64_t* device_basis,             //cantor basis (flattened, unused in simplified formula)
    int N,                              //current codeword size
    int basis_len                       //basis length (log2(N), unused in simplified formula)
) {
    int I = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(I >= N / 2) return;
    
    //indices for current[2I] and current[2I+1]
    int idx_even = 2 * I * FIELD_WORDS;
    int idx_odd = (2 * I + 1) * FIELD_WORDS;
    int idx_out = I * FIELD_WORDS;
    
    //load current[2I] and current[2I+1]
    uint64_t c_even[FIELD_WORDS], c_odd[FIELD_WORDS];
    for(int w = 0; w < FIELD_WORDS; w++) {
        c_even[w] = device_codeword[idx_even + w];
        c_odd[w] = device_codeword[idx_odd + w];
    }
    
    //compute difference: c_odd - c_even
    uint64_t diff[FIELD_WORDS];
    field_sub(diff, c_odd, c_even, FIELD_WORDS);
    
    //compute slope: alpha * diff (for cantor basis, denominator simplifies out)
    uint64_t slope[FIELD_WORDS];
    field_mul(slope, device_alpha, diff, FIELD_WORDS);
    
    //compute next[I] = c_even + slope
    uint64_t result[FIELD_WORDS];
    field_add(result, c_even, slope, FIELD_WORDS);
    
    //store result
    for(int w = 0; w < FIELD_WORDS; w++) {
        device_codeword_nxt[idx_out + w] = result[w];
    }
    
    //debugging - uncomment to write intermediate values
    //if(I == 0 && blockIdx.x == 0) {
    //    print_field_device("c_even", c_even);
    //    print_field_device("c_odd", c_odd);
    //    print_field_device("diff", diff);
    //    print_field_device("slope", slope);
    //    print_field_device("result", result);
    //}
}

//kernel to build merkle tree layer
//layer 0: just field elements (2^20)
//layers 1-15: field || hash (concat_words = 8)
//layers 16-20: just hash (hash_words = 4)
__global__ void compute_tree_layer(
    uint64_t* device_codeword_nxt,      //next folded codeword
    uint64_t* device_tree_layer,        //current tree layer
    uint64_t* device_tree_layer_nxt,    //next tree layer
    uint64_t* device_digest,            //hash outputs
    int N,                              //current size
    int layer_num                       //which layer (0-20)
) {
    int I = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(I >= N / 2) return;
    
    //determine word size based on layer
    int words_current = (layer_num == 0) ? FIELD_WORDS : 
                       (layer_num <= 15) ? CONCAT_WORDS : HASH_WORDS;
    int words_next = (layer_num < 15) ? CONCAT_WORDS : HASH_WORDS;
    
    int idx_left = 2 * I * words_current;
    int idx_right = (2 * I + 1) * words_current;
    int idx_out = I * words_next;
    int idx_hash = I * HASH_WORDS;
    
    //combine siblings for hashing
    uint64_t combined[2 * CONCAT_WORDS];  //max size
    int combined_size = 2 * words_current;
    
    for(int w = 0; w < words_current; w++) {
        combined[w] = device_tree_layer[idx_left + w];
        combined[words_current + w] = device_tree_layer[idx_right + w];
    }
    
    //hash the combined siblings using SHA3
    SHA3((uint8_t*)&device_digest[idx_hash], 
         (uint8_t*)combined, 
         combined_size * sizeof(uint64_t), 
         256);
    
    //build next layer element
    if(layer_num < 15) {
        //layers 1-14: store field || hash
        for(int w = 0; w < FIELD_WORDS; w++) {
            device_tree_layer_nxt[idx_out + w] = device_codeword_nxt[I * FIELD_WORDS + w];
        }
        for(int w = 0; w < HASH_WORDS; w++) {
            device_tree_layer_nxt[idx_out + FIELD_WORDS + w] = device_digest[idx_hash + w];
        }
    } else {
        //layers 15-20: store just hash
        for(int w = 0; w < HASH_WORDS; w++) {
            device_tree_layer_nxt[idx_out + w] = device_digest[idx_hash + w];
        }
    }
}

//host function to launch commit kernel (simplified for cantor basis)
void commit_launch(
    std::vector<FieldT>& codeword,           //current codeword
    std::vector<FieldT>& codeword_nxt,       //next codeword (output)
    FieldT& alpha,                           //challenge
    int N,                                   //current size
    int layer_num,                           //which layer (0-20)
    std::vector<std::vector<FieldT>>& merkle_tree_layers  //merkle tree storage
) {
    printf("commit_launch: N=%d, layer=%d\n", N, layer_num);
    
    int basis_len = (int)log2(N);
    
    //flatten codeword to uint64_t array
    uint64_t* flat_codeword = new uint64_t[N * FIELD_WORDS];
    for(int i = 0; i < N; i++) {
        //TODO: extract uint64_t words from FieldT element
        //this depends on libff::gf256 internal representation
        //for now, placeholder
        for(int w = 0; w < FIELD_WORDS; w++) {
            flat_codeword[i * FIELD_WORDS + w] = 0;  //replace with actual data
        }
    }
    
    //flatten alpha (only parameter we need)
    uint64_t flat_alpha[FIELD_WORDS] = {0};
    
    //allocate device memory (simplified - only need codeword and alpha)
    uint64_t *d_codeword, *d_codeword_nxt, *d_alpha;
    
    cudaMalloc(&d_codeword, N * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc(&d_codeword_nxt, (N/2) * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc(&d_alpha, FIELD_WORDS * sizeof(uint64_t));
    
    //copy to device
    cudaMemcpy(d_codeword, flat_codeword, N * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(d_alpha, flat_alpha, FIELD_WORDS * sizeof(uint64_t), cudaMemcpyHostToDevice);
    
    //launch kernel (simplified parameters)
    int threads = 256;
    int blocks = (N/2 + threads - 1) / threads;
    
    commit_kernel<<<blocks, threads>>>(
        d_codeword, d_codeword_nxt, d_alpha, 
        nullptr, nullptr,  // offset and basis unused
        N, basis_len
    );
    
    cudaDeviceSynchronize();
    
    //check for errors
    cudaError_t err = cudaGetLastError();
    if(err != cudaSuccess) {
        printf("kernel error: %s\n", cudaGetErrorString(err));
    }
    
    //copy back
    uint64_t* flat_codeword_nxt = new uint64_t[(N/2) * FIELD_WORDS];
    cudaMemcpy(flat_codeword_nxt, d_codeword_nxt, (N/2) * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    
    //unflatten to codeword_nxt
    codeword_nxt.resize(N/2);
    for(int i = 0; i < N/2; i++) {
        //TODO: construct FieldT from uint64_t words
        //codeword_nxt[i] = ...
    }
    
    //cleanup
    cudaFree(d_codeword);
    cudaFree(d_codeword_nxt);
    cudaFree(d_alpha);
    
    delete[] flat_codeword;
    delete[] flat_codeword_nxt;
    
    printf("commit_launch completed\n");
}
