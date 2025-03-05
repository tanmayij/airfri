// #ifndef COMMIT_KERNEL_H
// #define COMMIT_KERNEL_H

// #include <stdint.h>
// #include <stdio.h>
// #include "../include/field.cuh"
// #include "../include/hash.cuh"
// #define MAX_EVAL_BASIS_LEN 20
// #define MAX_FRI_PARAMETERS 16
// #define HASH_WORDS 4
// const size_t CONCAT_WORDS = FIELD_WORDS + HASH_WORDS;
// #ifndef HASH_SIZE
// #define HASH_SIZE (HASH_WORDS * sizeof(uint64_t))
// #endif
// //const size_t field_words = 4;
// __device__ void print_field_kernel(const char *label, const uint64_t *field, int field_words);
// __host__ void print_field_host(const char *label, const uint64_t *field, int field_words);

// __global__ void commit_kernel(
//     uint64_t *device_codeword, uint64_t *device_codeword_nxt,
//     uint64_t *device_alpha, uint64_t *device_offset,
//     uint64_t *device_denominator_inv, uint64_t *device_eval_basis,
//     uint64_t *device_temp1, uint64_t *device_temp2,
//     uint64_t *device_temp3, uint64_t *device_temp4,
//     uint64_t *device_temp5, uint64_t *device_alpha_offset, 
//     uint64_t *tree_initial_leaf, uint64_t *tree_concat_words,
//     int N, int basis_len
// );

// __global__ void merkle_kernel(
//     uint64_t *device_layer_hashes, 
//     uint64_t *device_merkle_root, 
//     uint64_t *tree_hash_words, 
//     int N
// );
// extern "C" {
//     void commit_launch(
//         uint64_t **codeword, uint64_t **codeword_nxt, 
//         uint64_t *alpha, uint64_t *offset, 
//         uint64_t denominator_inv, uint64_t *eval_basis, 
//         int N, uint64_t *root,
//         uint64_t *tree_initial_leaf, uint64_t *tree_concat_words, uint64_t *tree_hash_words
//     );
// }
// __host__ __device__ void i_th_ele_in_span(uint64_t *result, uint64_t *basis, int len_basis, int i);


// #endif // COMMIT_KERNEL_CUH


#ifndef COMMIT_KERNEL_H
#define COMMIT_KERNEL_H

#include <stdint.h>
#include <stdio.h>
#include "../include/field.cuh"
#include "../include/hash-host.cuh"
#include "../include/hash.cuh"
#define MAX_EVAL_BASIS_LEN 20
#define MAX_FRI_PARAMETERS 16
#define HASH_WORDS 4
const size_t CONCAT_WORDS = FIELD_WORDS + HASH_WORDS;
#ifndef HASH_SIZE
#define HASH_SIZE (HASH_WORDS * sizeof(uint64_t))
#endif
//const size_t field_words = 4;
__device__ void print_field_kernel(const char *label, const uint64_t *field, int field_words);
__host__ void print_field_host(const char *label, const uint64_t *field, int field_words);

__global__ void commit_kernel(
    uint64_t *device_codeword, uint64_t *device_codeword_nxt,
    uint64_t *device_alpha, uint64_t *device_offset,
    uint64_t *device_denominator_inv, uint64_t *device_eval_basis,
    uint64_t *device_temp1, uint64_t *device_temp2,
    uint64_t *device_temp3, uint64_t *device_temp4,
    uint64_t *device_temp5, uint64_t *device_alpha_offset, 
    uint64_t *device_layer_hashes, uint64_t *device_tree_layer, uint64_t *device_tree_layer_nxt, int N, int basis_len
);

__global__ void merkle_kernel(
    uint64_t *device_layer_hashes, 
    uint64_t *device_merkle_root, 
    uint64_t *device_tree_layer,
    uint64_t *device_tree_layer_nxt,
    uint64_t *device_combined_sibling_codewords,
    uint64_t *device_digest,
    uint64_t *device_combined_sibling_hashes,
    int N
);

__global__ void compute_tree_layers(uint64_t *device_codeword_nxt, uint64_t *device_layer_hashes, uint64_t *device_tree_layer,
    uint64_t *device_tree_layer_nxt, uint64_t *device_combined_sibling_codewords, uint64_t *device_concat_codeword_to_hash, uint64_t *device_digest, int N);

__global__ void compute_merkle_root_kernel(
    uint64_t *device_tree_layer,   
    uint64_t *device_merkle_root   
);

extern "C" {
    void commit_launch(
        uint64_t **codeword, uint64_t **codeword_nxt, 
        uint64_t *alpha, uint64_t *offset, 
        uint64_t denominator_inv, uint64_t *eval_basis, 
        int N, uint64_t *root, uint64_t **tree_layer, uint64_t **tree_layer_nxt, uint64_t ***tree, int last_round, bool is_last_round
    );

}
__host__ __device__ void i_th_ele_in_span(uint64_t *result, uint64_t *basis, int len_basis, int i);


#endif // COMMIT_KERNEL_CUH