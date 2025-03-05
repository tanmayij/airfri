//search the merkle tree for a requested value
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "../include/fri.cuh"
#define HASH_WORDS 4  

// Example: Define a sample tree with three layers (you will replace this with your actual tree)
// #define MAX_LAYERS 17
// #define MAX_NODES  131072
// #define CONCAT_WORDS 8

// uint64_t tree[MAX_LAYERS][MAX_NODES][CONCAT_WORDS];  

// Function to check if two hashes are equal
int hash_match(uint64_t *hash1, uint64_t *hash2) {
    for (int i = 0; i < HASH_WORDS; i++) {
        if (hash1[i] != hash2[i]) {
            return 0; // Mismatch found
        }
    }
    return 1; // Hashes are equal
}

// Function to search for a hash in the Merkle tree
void search_tree(uint64_t *target_hash) {
    for (int layer = 0; layer < 12; layer++) {
        for (int index = 0; index < 131072; index++) {
            if (hash_match(tree[layer][index], target_hash)) {
                printf("Hash found at Layer %d, Index %d\n", layer, index);
                return;
            }
        }
    }
    printf("Hash not found in tree\n");
}

int main() {
    // Sample input hash (replace with user input)
    uint64_t target_hash[HASH_WORDS] = {
        0x68b82b8ef32b7f3b, 0x4f80eccf06373c5d,
        0x14c92eaaa4a5ebf7, 0x078a635cdd2d6360
    };

    // Call the search function
    search_tree(target_hash);

    return 0;
}