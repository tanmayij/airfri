
#include <openssl/evp.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include "../include/merkle.cuh"
#include "../include/hash-host.cuh"
#include <stdio.h>
#define FIELD_WORDS 4
const size_t CONCAT_WORDS = 8;  // FIELD_WORDS + HASH_WORDS

void print_field_merkle(const char *label, const uint64_t *field, int field_words) {
    printf("%s: ", label);
    for (int i = 0; i < field_words; i++) {
        printf("%016llx ", field[i]);
    }
    printf("\n");
}

// Helper function to hash data using SHA3-256
void hash_sha3_256(const uint64_t *data, size_t len, uint64_t *out) {
    EVP_MD_CTX *mdctx = EVP_MD_CTX_new();
    const EVP_MD *md = EVP_sha3_256();

    EVP_DigestInit_ex(mdctx, md, NULL);
    EVP_DigestUpdate(mdctx, data, len);
    EVP_DigestFinal_ex(mdctx, (unsigned char *)out, NULL);

    EVP_MD_CTX_free(mdctx);
}

// Function to print the byte array
void print_bytes(const uint64_t *bytes, size_t size) {
    for (size_t i = 0; i < size; i++) {
        printf("%016llx", bytes[i]);
        if (i < size - 1) {
            printf(" ");
        }
    }
    printf("\n");
}

// Function to compute the Merkle root from leaf hashes
void merkle_commit(uint64_t **leaf_hashes, size_t num_leaves, uint64_t *out) {
    // Assert to ensure memory is allocated correctly and no corruption occurs
    for (size_t i = 0; i < num_leaves; i++) {
        assert(leaf_hashes[i] != NULL && "leaf_hashes[i] is not properly allocated");
        //printf("leaf_hashes[%zu]: %p\n", i, (void*)leaf_hashes[i]);  // Debug to print pointers of leaf_hashes
    }

    assert((num_leaves & (num_leaves - 1)) == 0); // num_leaves must be a power of two

    if (num_leaves == 1) {
        //printf("Computing Merkle root for single leaf\n");
        memcpy(out, leaf_hashes[0], HASH_SIZE);
        return;
    }

    size_t half = num_leaves / 2;
    uint64_t left_root[HASH_WORDS];
    uint64_t right_root[HASH_WORDS];

    //printf("Recursing left on %zu leaves\n", half);
    merkle_commit(leaf_hashes, half, left_root);
    //printf("Recursing right on %zu leaves\n", half);
    merkle_commit(leaf_hashes + half, half, right_root);

    uint64_t combined[2 * HASH_WORDS];
    memcpy(combined, left_root, HASH_SIZE);
    memcpy(combined + HASH_WORDS, right_root, HASH_SIZE);

    // Compute the hash of the concatenated hashes
    hash_sha3_256(combined, 2 * HASH_SIZE, out);
    //printf("Computed root for current level\n");
}

// Function to check if a point (hash or codeword) has already been sent
bool is_sent(uint64_t *point, uint64_t **sent_points, size_t sent_count) {
    for (size_t i = 0; i < sent_count; i++) {
        if (memcmp(point, sent_points[i], HASH_SIZE) == 0) {
            return true;
        }
    }
    return false;
}

//the below function computes the authentication path of the merkle tree for a given index. it "opens" the commitment (tree) at given indices. 
//first, the tree needs to be flattened. below function flattens the 3d tree to a 1d array where the first n elements (of FIELD_WORDS length) are the fist n
//codeword elements of the initial codeword. 
typedef struct {
    uint64_t ***tree;      
    size_t *row_sizes;     
    size_t *col_sizes;     
    size_t num_layers;     
} FRI_Tree;

typedef struct {
    uint64_t *sibling_node;
    size_t sibling_node_size;
} AuthPathElement;

void get_fri_auth_path(FRI_Tree *fri_tree, size_t row_size, size_t col_size, AuthPathElement **auth_path, size_t *path_length) {
    size_t current_row = leaf_row;
    size_t current_col = leaf_col;
    size_t num_layers = fri_tree->num_layers;

    *auth_path = malloc(num_layers * sizeof(AuthPathElement));
    *path_length = 0;

    for (size_t layer = 0; layer < num_layers - 1; layer++) { // Traverse up to the root
        size_t row_size = fri_tree->row_sizes[layer];
        size_t col_size = fri_tree->col_sizes[layer];

        size_t sibling_row = (current_row % 2 == 0) ? current_row + 1 : current_row - 1;
        size_t sibling_col = (current_col % 2 == 0) ? current_col + 1 : current_col - 1;

        sibling_row = (sibling_row < row_size) ? sibling_row : current_row;
        sibling_col = (sibling_col < col_size) ? sibling_col : current_col;

        (*auth_path)[*path_length].sibling_data = fri_tree->tree[layer][sibling_row][sibling_col];
        (*auth_path)[*path_length].size = fri_tree->col_sizes[layer]; // Assume full column size here
        (*path_length)++;
        current_row /= 2;
        current_col /= 2;
    }
} 

uint64_t *flatten_Merkletree(uint64_t ***tree, uint64_t *output_flattened_tree){
    
}
int merkle_open(uint64_t )