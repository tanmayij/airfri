// #ifndef MERKLE_H
// #define MERKLE_H

// #include <openssl/evp.h>
// #include <stdint.h>
// #include <stdlib.h>
// #include <string.h>
// #include <assert.h>
// #include "../include/hash-host.cuh"
// #include <stdio.h>

// //size_t field_words = 4; //for 192
// #define HASH_WORDS 4 //64-bit words
// #define HASH_SIZE (HASH_WORDS * sizeof(uint64_t))

// //Helper function to hash data using SHA3-256
// void hash_sha3_256(const uint64_t *data, size_t len, uint64_t *out);

// //Function to print the byte array
// void print_bytes(const uint64_t *bytes, size_t size);

// //Function to compute the Merkle root from leaf hashes
// void merkle_commit(uint64_t **leaf_hashes, size_t num_leaves, uint64_t *out);


// void merkle_open(
//     uint64_t ***codewords,     // 3D array for codewords at each layer: codewords[layer][index]
//     size_t num_leaves,         // Initial number of leaves (size of the first codeword)
//     size_t leaf_index,              // Index of the leaf to be proven
//     uint64_t **proof_path,     // Output array for the proof path (each entry is a hash-sized block)
//     size_t *proof_len,         // Output length of the proof path
//     size_t field_words         // Number of 64-bit words in each field element
// );

// bool merkle_verify(
//     uint64_t *root_hash,          // Known root hash of the Merkle tree
//     uint64_t *leaf,               // Initial codeword element (leaf node)
//     size_t leaf_idx,              // Index of the leaf element
//     uint64_t **proof_path,        // Proof path to verify
//     size_t proof_len,             // Length of the proof path
//     size_t field_words            // Number of words in each field element
// );

// #endif // MERKLE_H

#ifndef MERKLE_H
#define MERKLE_H

#include <openssl/evp.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../include/hash-host.cuh"
#include <stdio.h>

//size_t field_words = 4; //for 192
#define HASH_WORDS 4 //64-bit words
#define HASH_SIZE (HASH_WORDS * sizeof(uint64_t))

//Helper function to hash data using SHA3-256
void hash_sha3_256(const uint64_t *data, size_t len, uint64_t *out);

//Function to print the byte array
void print_bytes(const uint64_t *bytes, size_t size);

//Function to compute the Merkle root from leaf hashes
void merkle_commit(uint64_t **leaf_hashes, size_t num_leaves, uint64_t *out);


void merkle_open(uint64_t **auth_path, int leaf_idx, size_t *proof_len, uint64_t *host_concatenated_tree);

int merkle_verify(uint64_t *root, size_t index, uint64_t **auth_path, size_t proof_len, uint64_t *leaf);

#endif // MERKLE_H