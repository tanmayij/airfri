
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

void merkle_open(uint64_t **auth_path, int leaf_idx, size_t *proof_len, uint64_t ***tree, int layer) {
    int current_index = leaf_idx;
    *proof_len = 0;

    for (int i = layer; i < 16; i++) {
        int sibling_index = (current_index % 2 == 0) ? current_index + 1 : current_index - 1;
        size_t sibling_size;

        if (i == 0) {
            // For layer 0, only copy the sibling (FIELD_WORDS size)
            sibling_size = FIELD_WORDS;
            auth_path[*proof_len] = (uint64_t *)malloc(sibling_size * sizeof(uint64_t));
            if (!auth_path[*proof_len]) {
                fprintf(stderr, "Memory allocation failed for auth_path at proof_len %zu\n", *proof_len);
                exit(1);
            }
            memcpy(auth_path[*proof_len], tree[i][sibling_index], sibling_size * sizeof(uint64_t));
            printf("trying to set a breakpoint here tanjan.\n");

        } else if (i > 0 && i < 11) {
            // For layers 1 to 12, include the codeword element (FIELD_WORDS) and sibling
            sibling_size = FIELD_WORDS + HASH_WORDS;
            auth_path[*proof_len] = (uint64_t *)malloc((FIELD_WORDS + sibling_size) * sizeof(uint64_t));
            if (!auth_path[*proof_len]) {
                fprintf(stderr, "Memory allocation failed for auth_path at proof_len %zu\n", *proof_len);
                exit(1);
            }
            // Copy the codeword element from the current index (skipping the hash part)
            memcpy(auth_path[*proof_len], tree[i][current_index], FIELD_WORDS * sizeof(uint64_t));
            // Copy the sibling element in full
            memcpy(auth_path[*proof_len] + FIELD_WORDS, tree[i][sibling_index], sibling_size * sizeof(uint64_t));

        } else {
            // For layers 13 to 16, only include the sibling (HASH_WORDS size)
            sibling_size = HASH_WORDS;
            auth_path[*proof_len] = (uint64_t *)malloc(sibling_size * sizeof(uint64_t));
            if (!auth_path[*proof_len]) {
                fprintf(stderr, "Memory allocation failed for auth_path at proof_len %zu\n", *proof_len);
                exit(1);
            }
            memcpy(auth_path[*proof_len], tree[i][sibling_index], sibling_size * sizeof(uint64_t));
        }

        (*proof_len)++;
        current_index /= 2;  // Move up the tree for the next layer
        //printf("current index is: %d\n", current_index);
    }
}

// int merkle_verify(
//     uint64_t *root,
//     size_t leaf_idx,
//     uint64_t **auth_path,
//     size_t proof_len,
//     uint64_t *leaf,
//     int layer
// ) {
//     uint64_t current_hash[HASH_WORDS];
//     uint64_t combined[2 * FIELD_WORDS + HASH_WORDS];
    
//     // Initialize with the leaf value
//     memcpy(current_hash, leaf, FIELD_WORDS * sizeof(uint64_t));

//     for (int i = 0; i < proof_len; ++i) {
//         // Layer 0: only FIELD_WORDS in auth_path
//         if (i == 0) {
//             // Combine leaf with auth_path[0] for hash calculation
//             memcpy(combined, current_hash, FIELD_WORDS * sizeof(uint64_t));
//             memcpy(combined + FIELD_WORDS, auth_path[i], FIELD_WORDS * sizeof(uint64_t));
//             SHA3_host((uint8_t *)current_hash, (uint8_t *)combined, 2 * FIELD_WORDS * sizeof(uint64_t), 256);
//         }
//         // Layers 1-12: auth_path contains FIELD_WORDS + HASH_WORDS + FIELD_WORDS
//         else if (i > 0 && i < 12) {
//             // Extract first FIELD_WORDS for combination
//             memcpy(combined, auth_path[i], FIELD_WORDS * sizeof(uint64_t));
//             // Concat previous hash to this first part
//             memcpy(combined + FIELD_WORDS, current_hash, HASH_WORDS * sizeof(uint64_t));
//             // Concat remaining part of auth_path (HASH_WORDS + FIELD_WORDS)
//             memcpy(combined + FIELD_WORDS + HASH_WORDS, auth_path[i] + FIELD_WORDS, (HASH_WORDS + FIELD_WORDS) * sizeof(uint64_t));
//             SHA3_host((uint8_t *)current_hash, (uint8_t *)combined, (2 * FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t), 256);
//         }
//         // Layers 13 and above: auth_path contains only HASH_WORDS
//         else {
//             // Concat current hash with auth_path[i] for the final hash
//             memcpy(combined, current_hash, HASH_WORDS * sizeof(uint64_t));
//             memcpy(combined + HASH_WORDS, auth_path[i], HASH_WORDS * sizeof(uint64_t));
//             SHA3_host((uint8_t *)current_hash, (uint8_t *)combined, 2 * HASH_WORDS * sizeof(uint64_t), 256);
//         }
//     }
//     printf("Computed Merkle Root: ");
//     for (int i = 0; i < HASH_WORDS; i++) {
//         printf("%016lx ", current_hash[i]);
//     }
//     printf("\n");
//     // Final comparison with the provided root
//     return memcmp(root, root, HASH_WORDS * sizeof(uint64_t)) == 0;
// }

// int merkle_verify(
//     uint64_t *root,
//     size_t leaf_idx,
//     uint64_t **auth_path,
//     size_t proof_len,
//     uint64_t *leaf,
//     int layer
// ) { 
//     uint64_t hash[HASH_WORDS];
//     if(layer == 0){
//         uint64_t current_element[FIELD_WORDS];
//         memcpy(current_element, leaf, FIELD_WORDS * sizeof(uint64_t));

//         int sibling_size = FIELD_WORDS;
//         uint64_t combined[FIELD_WORDS + sibling_size];
//         memcpy(combined, current_element, FIELD_WORDS * sizeof(uint64_t));
//         memcpy(combined + FIELD_WORDS, auth_path_a[0], FIELD_WORDS * sizeof(uint64_t));
//         SHA3_host((uint8_t *)hash, (uint8_t *)combined, (FIELD_WORDS + sibling_size) * sizeof(uint64_t));
//     }
// }

int merkle_verify(
    uint64_t *root,           // Expected Merkle root to verify against
    size_t leaf_idx,          // Index of the leaf being verified
    uint64_t **auth_path,     // Authentication path (sibling hashes) for the leaf
    size_t proof_len,         // Length of the proof path
    uint64_t *leaf,            // Initial leaf (codeword element) to start the verification
    int layer
) {
    uint64_t current_hash[HASH_WORDS];
    memcpy(current_hash, leaf, FIELD_WORDS * sizeof(uint64_t));

    for (size_t i = layer; i < proof_len; i++) {
        uint64_t combined[2 * (FIELD_WORDS + HASH_WORDS)]; // Buffer for concatenated data
        size_t combined_size;

        if (i == 0) {
            // Layer 0: Concatenate the codeword element with sibling, both FIELD_WORDS size
            combined_size = 2 * FIELD_WORDS;
            memcpy(combined, current_hash, FIELD_WORDS * sizeof(uint64_t));
            memcpy(combined + FIELD_WORDS, auth_path[i], FIELD_WORDS * sizeof(uint64_t));
        } else if (i < 11) {
            // Layers 1-12: Concatenate current hash with the first part of auth_path and the full sibling
            combined_size = FIELD_WORDS + (FIELD_WORDS + HASH_WORDS);
            memcpy(combined, current_hash, HASH_WORDS * sizeof(uint64_t));
            memcpy(combined + HASH_WORDS, auth_path[i], (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
        } else {
            // Layers 13-16: Only hash concatenation with sibling, both HASH_WORDS size
            combined_size = 2 * HASH_WORDS;
            memcpy(combined, current_hash, HASH_WORDS * sizeof(uint64_t));
            memcpy(combined + HASH_WORDS, auth_path[i], HASH_WORDS * sizeof(uint64_t));
        }

        // Hash the concatenated result to get the new current hash
        SHA3_host((uint8_t *)current_hash, (uint8_t *)combined, combined_size * sizeof(uint64_t), 256);
    }
    printf("Computed Merkle Root: ");
            for (int i = 0; i < HASH_WORDS; i++) {
                printf("%016lx ", current_hash[i]);
            }
            printf("\n");
    // Compare the computed root with the expected root
    return memcmp(root, current_hash, HASH_WORDS * sizeof(uint64_t)) == 0;
}