#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <assert.h>
#include <string.h>
#include <cuda_runtime.h>
#include "../include/commit-launch-merkle.cuh"
#include "../include/merkle.cuh"
#include "../include/field.cuh"

// Define constants for the test
#define FIELD_WORDS 4
#define HASH_WORDS 4

// Initialize a fixed eval_basis and offset for testing consistency
uint64_t eval_basis[] = {0x1, 0x2, 0x4, 0x8};  // Sample values
uint64_t offset = 0x1;

// Function to print field elements for debugging
void print_field(const char *label, const uint64_t *field, size_t field_words) {
    printf("%s: ", label);
    for (size_t i = 0; i < field_words; i++) {
        printf("%016lx ", field[i]);
    }
    printf("\n");
}

// Function to run the Merkle Tree commitment and verification process
void test_merkle_tree_commitment_and_verification() {
    // Step 1: Initialize an initial codeword with 8 elements
    size_t N = 8;
    uint64_t ***codewords = (uint64_t ***)malloc(3 * sizeof(uint64_t **));  // Array to store all codewords

    // Allocate memory for each round's codewords array and initialize the first round
    for (int round = 0; round < 3; round++) {
        size_t round_size = N >> round;  // Size halves each round
        codewords[round] = (uint64_t **)malloc(round_size * sizeof(uint64_t *));
        
        if (round == 0) {
            for (size_t i = 0; i < round_size; i++) {
                codewords[round][i] = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
                for (size_t j = 0; j < FIELD_WORDS; j++) {
                    codewords[round][i][j] = i + j;  // Simplified assignment for testing
                }
            }
        }
    }

    // Array to store the Merkle root after each round
    uint64_t *root = (uint64_t *)malloc(HASH_WORDS * sizeof(uint64_t));
    uint64_t *tree = (uint64_t *)malloc(4 * (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));

    // Step 2: Run commit_launch for 3 rounds, each time halving the codeword size
    for (int round = 0; round < 3; round++) {
        printf("\n--- Round %d ---\n", round + 1);
        size_t new_N = N / 2;

        uint64_t **codeword_nxt = (uint64_t **)malloc(new_N * sizeof(uint64_t *));
        // for (size_t i = 0; i < new_N; i++) {
        //     codeword_nxt[i] = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
        // }

        uint64_t alpha[FIELD_WORDS] = {0x2};  // Example alpha value
        uint64_t denominator_inv = 0x1;  // Example denominator inverse
        
        commit_launch(codewords[round], codeword_nxt, alpha, &offset, denominator_inv, eval_basis, N, root, tree);

        // Store the new codeword_nxt into codewords for the next round
        if (round < 2) {
            for (size_t i = 0; i < new_N; i++) {
                codewords[round + 1][i] = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
                memcpy(codewords[round + 1][i], codeword_nxt[i], FIELD_WORDS * sizeof(uint64_t));
            }
        }

        // Print Merkle root for each round
        printf("Merkle Root after round %d: ", round + 1);
        print_field("Root", root, HASH_WORDS);

        // Free codeword_nxt as it has been copied to codewords[round + 1]
        for (size_t i = 0; i < new_N; i++) {
            free(codeword_nxt[i]);
        }
        free(codeword_nxt);
        N = new_N;
    }

    // Step 3: Verify the Merkle root and authentication paths for the final codeword
    printf("\n--- Final Merkle Tree Verification ---\n");
    for (size_t i = 0; i < N; i++) {
        int proof_len = 0;
        uint64_t *auth_path = (uint64_t *)malloc(4 * CONCAT_WORDS * sizeof(uint64_t));
        for (size_t j = 0; j < N; j++) {
            auth_path[j] = (uint64_t)malloc(HASH_SIZE);
        }

        uint64_t *leaf = *codewords[2];
        for(int k =0; k < FIELD_WORDS; k++){
            printf("%016lx ", leaf[k]);
        }
        merkle_open(auth_path, i, proof_len, tree);

        int result = merkle_verify(root, i, auth_path, proof_len, leaf);
        printf("Verification for leaf %zu: %s\n", i, result ? "Success" : "Failure");

        
        free(auth_path);
    }

    // Cleanup all codewords arrays
    for (int round = 0; round < 3; round++) {
        size_t round_size = 8 >> round;
        for (size_t i = 0; i < round_size; i++) {
            free(codewords[round][i]);
        }
        free(codewords[round]);
    }
    free(codewords);

    printf("\nMerkle Tree Commitment and Verification Test Completed.\n");
}

int main() {
    test_merkle_tree_commitment_and_verification();
    return 0;
}