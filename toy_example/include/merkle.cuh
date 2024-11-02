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


void merkle_open(uint64_t *proof_path, int leaf_idx, int proof_len, uint64_t* host_concatenated_tree);

int merkle_verify(uint64_t *root, size_t index, uint64_t *proof_path, size_t proof_len, uint64_t *leaf);

#endif // MERKLE_H