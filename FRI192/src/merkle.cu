
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

void merkle_open(uint64_t **auth_path, int leaf_idx, size_t *proof_len, uint64_t ***tree) {
    int current_index = leaf_idx;
    *proof_len = 0;
    int layer = 0;

    while (layer < 17) {
        int sibling_index = (current_index % 2 == 0) ? current_index + 1 : current_index - 1;
        size_t sibling_size;
        if (layer == 0) {
            sibling_size = FIELD_WORDS; 
        } else if (layer < 13) {
            sibling_size = FIELD_WORDS + HASH_WORDS; 
        } else {
            sibling_size = HASH_WORDS; 
        }
        auth_path[*proof_len] = (uint64_t *)malloc((HASH_WORDS + sibling_size) * sizeof(uint64_t));
        if (!auth_path[*proof_len]) {
            fprintf(stderr, "Memory allocation failed for auth_path at proof_len %zu\n", *proof_len);
            exit(1);
        }
        memcpy(auth_path[*proof_len], tree[layer][current_index], HASH_WORDS * sizeof(uint64_t)); //pick the hash from each index
        memcpy(auth_path[*proof_len] + HASH_WORDS, tree[layer][sibling_index], sibling_size * sizeof(uint64_t));      //add sibling to it
        
        (*proof_len)++;
        current_index /= 2;
        layer++;
    }
}

// int merkle_verify(
//     uint64_t *root,           // Expected Merkle root to verify against
//     size_t leaf_idx,          // Index of the leaf being verified
//     uint64_t **auth_path,     // Authentication path (sibling hashes) for the leaf
//     size_t proof_len,         // Length of the proof path
//     uint64_t *leaf            // Initial leaf (codeword element) to start the verification
// ) {
//     uint64_t computed_hash[HASH_WORDS];
//     memcpy(computed_hash, leaf, FIELD_WORDS * sizeof(uint64_t)); 

//     size_t current_idx = leaf_idx;

//     for (size_t i = 0; i < proof_len; ++i) {
//         uint64_t combined[FIELD_WORDS + 2 * HASH_WORDS];
//         size_t combined_size;

//         if (i == 0) {
//             combined_size = FIELD_WORDS + HASH_WORDS;
//             if (current_idx % 2 == 0) {
//                 memcpy(combined, computed_hash, FIELD_WORDS * sizeof(uint64_t));
//                 memcpy(combined + FIELD_WORDS, auth_path[i], HASH_WORDS * sizeof(uint64_t));
//             } else {
//                 memcpy(combined, auth_path[i], HASH_WORDS * sizeof(uint64_t));
//                 memcpy(combined + HASH_WORDS, computed_hash, FIELD_WORDS * sizeof(uint64_t));
//             }
//         } else if (i < 13) {
//             combined_size = FIELD_WORDS + HASH_WORDS;
//             if (current_idx % 2 == 0) {
//                 memcpy(combined, computed_hash, HASH_WORDS * sizeof(uint64_t));
//                 memcpy(combined + HASH_WORDS, auth_path[i], (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
//             } else {
//                 memcpy(combined, auth_path[i], (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
//                 memcpy(combined + FIELD_WORDS + HASH_WORDS, computed_hash, HASH_WORDS * sizeof(uint64_t));
//             }
//         } else {
//             combined_size = 2 * HASH_WORDS;
//             if (current_idx % 2 == 0) {
//                 memcpy(combined, computed_hash, HASH_WORDS * sizeof(uint64_t));
//                 memcpy(combined + HASH_WORDS, auth_path[i], HASH_WORDS * sizeof(uint64_t));
//             } else {
//                 memcpy(combined, auth_path[i], HASH_WORDS * sizeof(uint64_t));
//                 memcpy(combined + HASH_WORDS, computed_hash, HASH_WORDS * sizeof(uint64_t));
//             }
//         }
//         SHA3_host((uint8_t *)computed_hash, (uint8_t *)combined, combined_size * sizeof(uint64_t), 256);
//         current_idx /= 2;
//     }
//     return memcmp(computed_hash, root, HASH_WORDS * sizeof(uint64_t)) == 0;
// }

int merkle_verify(
    uint64_t *root,           // Expected Merkle root to verify against
    size_t leaf_idx,          // Index of the leaf being verified
    uint64_t **auth_path,     // Authentication path (sibling hashes) for the leaf
    size_t proof_len,         // Length of the proof path
    uint64_t *leaf            // Initial leaf (codeword element) to start the verification
) {
    uint64_t current_codeword[FIELD_WORDS];
    memcpy(current_codeword, leaf, FIELD_WORDS * sizeof(uint64_t)); 

    size_t current_idx = leaf_idx;

    for (size_t i = 0; i < proof_len; ++i) {
        size_t combined_size;

        if (i == 0) {
            uint64_t combined[2 * FIELD_WORDS];  // Buffer for concatenation
            combined_size = 2 * FIELD_WORDS;
            if (current_idx % 2 == 0) {
                // Leaf is on the left
                memcpy(combined, auth_path[i], 2 * FIELD_WORDS * sizeof(uint64_t));
            } else {
                // Leaf is on the right
                memcpy(combined, auth_path[i], FIELD_WORDS * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS, auth_path[i + FIELD_WORDS], FIELD_WORDS * sizeof(uint64_t));
            }
        } else if (i < 13) {
            // Intermediate layers with FIELD_WORDS + HASH_WORDS
            combined_size = FIELD_WORDS + HASH_WORDS;
            if (current_idx % 2 == 0) {
                // Computed hash on the left, sibling on the right
                memcpy(combined, computed_hash, HASH_WORDS * sizeof(uint64_t));
                memcpy(combined + HASH_WORDS, auth_path[i], (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
            } else {
                // Sibling on the left, computed hash on the right
                memcpy(combined, auth_path[i], (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS + HASH_WORDS, computed_hash, HASH_WORDS * sizeof(uint64_t));
            }
        } else {
            // Final layers with only HASH_WORDS
            combined_size = 2 * HASH_WORDS;
            if (current_idx % 2 == 0) {
                // Computed hash on the left, sibling on the right
                memcpy(combined, computed_hash, HASH_WORDS * sizeof(uint64_t));
                memcpy(combined + HASH_WORDS, auth_path[i], HASH_WORDS * sizeof(uint64_t));
            } else {
                // Sibling on the left, computed hash on the right
                memcpy(combined, auth_path[i], HASH_WORDS * sizeof(uint64_t));
                memcpy(combined + HASH_WORDS, computed_hash, HASH_WORDS * sizeof(uint64_t));
            }
        }

        // Hash the concatenated result
        SHA3_host((uint8_t *)computed_hash, (uint8_t *)combined, combined_size * sizeof(uint64_t), 256);

        // Move up to the next layer
        current_idx /= 2;
    }
    printf("Computed Merkle Root: ");
    for (int i = 0; i < HASH_WORDS; i++) {
        printf("%016lx ", computed_hash[i]);
    }
    printf("\n");
    // Check if computed hash matches the Merkle root
    return memcmp(computed_hash, root, HASH_WORDS * sizeof(uint64_t)) == 0;
}