
#include <openssl/evp.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <assert.h>
#include "../include/merkle.cuh"
#include "../include/hash-host.cuh"
#include <stdio.h>
#define FIELD_WORDS 5
const size_t CONCAT_WORDS = 9;  // FIELD_WORDS + HASH_WORDS

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

// void merkle_open(
//     merkleTree *tree,       // Complete Merkle tree array
//     int leaf_idx,           // Index of the leaf to open
//     int N,                  // Initial number of codewords
//     int current_layer,       // Target layer up to which the path is constructed (usually 16 for the root)
//     uint64_t *auth_path,    // Array to store the authentication path
//     size_t *auth_path_len   // Length of auth_path (used for tracking position)
// ) {
//     int current_idx = leaf_idx; // Start from the leaf node
//     int target_layer = 17;    
//     *auth_path_len = 0;          // Initialize the path length

//     while (current_layer < target_layer) {
//         int sibling_idx = (current_idx % 2 == 0) ? current_idx + 1 : current_idx - 1;
        
//         if (current_layer == 0) {
//             // Layer 0: Copy the sibling codeword only
//             memcpy(&auth_path[*auth_path_len], tree[sibling_idx].element, FIELD_WORDS * sizeof(uint64_t));
//             *auth_path_len += FIELD_WORDS;
//         } else if (current_layer <= 11) {
//             // Layers 1 to 11: Copy current index's codeword, sibling's codeword, and sibling's hash
//             memcpy(&auth_path[*auth_path_len], tree[current_idx].element, FIELD_WORDS * sizeof(uint64_t));
//             *auth_path_len += FIELD_WORDS;

//             memcpy(&auth_path[*auth_path_len], tree[sibling_idx].element, FIELD_WORDS * sizeof(uint64_t));
//             *auth_path_len += FIELD_WORDS;

//             memcpy(&auth_path[*auth_path_len], tree[sibling_idx].hash, HASH_WORDS * sizeof(uint64_t));
//             *auth_path_len += HASH_WORDS;
//         } else {
//             // Layers > 11: Only the sibling's hash
//             memcpy(&auth_path[*auth_path_len], tree[sibling_idx].hash, HASH_WORDS * sizeof(uint64_t));
//             *auth_path_len += HASH_WORDS;
//         }

//         // Move to the parent for the next loop iteration
//         current_idx /= 2;
//         current_layer++;
//     }
// }
void merkle_open(
    int leaf_idx,                 // Index of the leaf node for which we need the authentication path
    uint64_t **codewords,         // 2D array of codewords at each layer
    uint64_t *auth_path,          // Output array for storing the authentication path
    size_t *auth_path_len,        // Pointer to track the length of auth_path
    int N                          // Total number of leaves at the current layer
) {
    size_t I = leaf_idx;
    int layer = 0;
    *auth_path_len = 0;  // Initialize auth path length

    uint64_t combined[CONCAT_WORDS];  // Array to hold concatenated codeword and hash
    uint8_t digest[HASH_SIZE];        // For SHA3-256 hash results

    // Traverse layers until reaching the root
    while (N > 1) {
        int sibling_idx = (I % 2 == 0) ? I + 1 : I - 1;
        int idx1 = I * FIELD_WORDS;
        int idx2 = sibling_idx * FIELD_WORDS;

        if (layer == 0) {
            // Layer 0: Push only the sibling codeword element
            memcpy(&auth_path[*auth_path_len], &codewords[layer][idx2], FIELD_WORDS * sizeof(uint64_t));
            *auth_path_len += FIELD_WORDS;
        } 
        else if (layer > 0 && layer <= 11) {
            // Layers 1 to 11: Push current codeword, and concatenate sibling's hash || codeword element
            memcpy(&auth_path[*auth_path_len], &codewords[layer][idx1], FIELD_WORDS * sizeof(uint64_t));
            *auth_path_len += FIELD_WORDS;

            // Concatenate sibling's hash and codeword, then push without hashing
            memcpy(combined, &codewords[layer][idx2], FIELD_WORDS * sizeof(uint64_t));
            SHA3_host(digest, (uint8_t *)combined, FIELD_WORDS * sizeof(uint64_t), 256);
            memcpy(combined + FIELD_WORDS, digest, HASH_WORDS * sizeof(uint64_t));
            memcpy(&auth_path[*auth_path_len], combined, CONCAT_WORDS * sizeof(uint64_t));
            *auth_path_len += CONCAT_WORDS;
        } 
        else {
            // Layers > 11: Only the sibling hash
            SHA3_host(digest, (uint8_t *)&codewords[layer][idx2], FIELD_WORDS * sizeof(uint64_t), 256);
            memcpy(&auth_path[*auth_path_len], digest, HASH_WORDS * sizeof(uint64_t));
            *auth_path_len += HASH_WORDS;
        }

        // Calculate the hash for the current layer to be used in the next iteration
        if (layer < 12) {
            memcpy(combined, &codewords[layer][idx1], FIELD_WORDS * sizeof(uint64_t));
            memcpy(combined + FIELD_WORDS, &codewords[layer][idx2], FIELD_WORDS * sizeof(uint64_t));
            SHA3_host(digest, (uint8_t *)combined, (FIELD_WORDS * 2) * sizeof(uint64_t), 256);
            memcpy(&codewords[layer + 1][I], digest, HASH_WORDS * sizeof(uint64_t));
        }

        // Move to the parent layer
        I /= 2;
        N /= 2;
        layer++;
    }
}
// void merkle_open(uint64_t **auth_path, int leaf_idx, size_t *proof_len, uint64_t *host_concatenated_tree) {
//     int current_index = leaf_idx;
//     *proof_len = 0;  // Initialize proof_len to 0 and increment as we add layers

//     // Traverse up the tree and collect sibling hashes
//     while (current_index > 0) {
//         int sibling_index;

//         // Determine sibling index: if current index is even, sibling is current_index + 1; if odd, it's current_index - 1
//         if (current_index % 2 == 0) {
//             sibling_index = current_index + 1;
//         } else {
//             sibling_index = current_index - 1;
//         }

//         // Allocate space for the sibling hash in auth_path
//         auth_path[*proof_len] = (uint64_t *)malloc(CONCAT_WORDS * sizeof(uint64_t));

//         // Copy the sibling hash from host_concatenated_tree to auth_path for the current layer
//         memcpy(auth_path[*proof_len], 
//                host_concatenated_tree + sibling_index * CONCAT_WORDS, 
//                CONCAT_WORDS * sizeof(uint64_t));

//         // Increment proof_len since we're adding one layer to the proof path
//         (*proof_len)++;

//         // Move up one level in the tree (integer division to find parent index)
//         current_index /= 2;
//     }
// }

// bool merkle_verify(
//     const uint64_t *leaf_element, // The element at the leaf
//     int leaf_idx,                 // Index of the leaf in the tree
//     const uint64_t *root_hash,    // Expected root hash for verification
//     const uint64_t *auth_path,    // Authentication path
//     size_t auth_path_len          // Length of the authentication path
// ) {
//     uint64_t computed_hash[HASH_WORDS];        // Buffer for the computed hash at each level
//     uint64_t concat_buffer[2 * FIELD_WORDS];   // Buffer for concatenating elements and hashes
//     int current_idx = leaf_idx;                // Start from the leaf index
//     int current_layer = 0;                     // Layer counter in the Merkle tree
//     size_t path_pos = 0;                       // Position in the auth_path array

//     // Traverse up the Merkle tree layers, applying verification rules per layer
//     while (current_layer <= 16 && path_pos < auth_path_len) {
//         uint64_t sibling_hash[HASH_WORDS];
//         uint64_t sibling_element[FIELD_WORDS];

//         // Handle concatenation and hashing per layer
//         if (current_layer == 0) {
//             // Layer 0: Concatenate leaf_element with sibling codeword from auth_path
//             memcpy(sibling_element, &auth_path[path_pos], FIELD_WORDS * sizeof(uint64_t));
//             path_pos += FIELD_WORDS;

//             // Concatenate and hash to compute the parent node
//             if (current_idx % 2 == 0) {
//                 memcpy(concat_buffer, leaf_element, FIELD_WORDS * sizeof(uint64_t));
//                 memcpy(concat_buffer + FIELD_WORDS, sibling_element, FIELD_WORDS * sizeof(uint64_t));
//             } else {
//                 memcpy(concat_buffer, sibling_element, FIELD_WORDS * sizeof(uint64_t));
//                 memcpy(concat_buffer + FIELD_WORDS, leaf_element, FIELD_WORDS * sizeof(uint64_t));
//             }
//             SHA3_host((uint8_t *)computed_hash, (uint8_t *)concat_buffer, 2 * FIELD_WORDS * sizeof(uint64_t), 256);

//         } else if (current_layer <= 11) {
//             // Layers 1 to 11: Process current codeword, sibling codeword, and sibling hash
//             uint64_t current_codeword[FIELD_WORDS];
//             memcpy(current_codeword, &auth_path[path_pos], FIELD_WORDS * sizeof(uint64_t));
//             path_pos += FIELD_WORDS;

//             memcpy(sibling_element, &auth_path[path_pos], FIELD_WORDS * sizeof(uint64_t));
//             path_pos += FIELD_WORDS;

//             memcpy(sibling_hash, &auth_path[path_pos], HASH_WORDS * sizeof(uint64_t));
//             path_pos += HASH_WORDS;

//             // Concatenate values depending on node position, then hash to compute parent node
//             if (current_idx % 2 == 0) {
//                 memcpy(concat_buffer, current_codeword, FIELD_WORDS * sizeof(uint64_t));
//                 memcpy(concat_buffer + FIELD_WORDS, sibling_hash, HASH_WORDS * sizeof(uint64_t));
//             } else {
//                 memcpy(concat_buffer, sibling_element, FIELD_WORDS * sizeof(uint64_t));
//                 memcpy(concat_buffer + FIELD_WORDS, computed_hash, HASH_WORDS * sizeof(uint64_t));
//             }
//             SHA3_host((uint8_t *)computed_hash, (uint8_t *)concat_buffer, (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t), 256);

//         } else {
//             // Layers > 11: Only sibling hashes
//             memcpy(sibling_hash, &auth_path[path_pos], HASH_WORDS * sizeof(uint64_t));
//             path_pos += HASH_WORDS;

//             // Concatenate sibling hash only, then hash to compute parent node
//             if (current_idx % 2 == 0) {
//                 memcpy(concat_buffer, computed_hash, HASH_WORDS * sizeof(uint64_t));
//                 memcpy(concat_buffer + HASH_WORDS, sibling_hash, HASH_WORDS * sizeof(uint64_t));
//             } else {
//                 memcpy(concat_buffer, sibling_hash, HASH_WORDS * sizeof(uint64_t));
//                 memcpy(concat_buffer + HASH_WORDS, computed_hash, HASH_WORDS * sizeof(uint64_t));
//             }
//             SHA3_host((uint8_t *)computed_hash, (uint8_t *)concat_buffer, 2 * HASH_WORDS * sizeof(uint64_t), 256);
//         }

//         // Move to the parent for the next iteration
//         current_idx /= 2;
//         current_layer++;
//     }

//     // Compare the computed root with the expected root hash
//     return (memcmp(computed_hash, root_hash, HASH_WORDS * sizeof(uint64_t)) == 0);
// }
bool merkle_verify(
    uint64_t *leaf_codeword,     // Codeword element at the specified leaf
    int leaf_idx,                // Index of the leaf to verify
    uint64_t *auth_path,         // Authentication path generated by merkle_open
    size_t auth_path_len,        // Length of the authentication path
    uint64_t *expected_root      // Expected Merkle root for verification
) {
    uint64_t current_hash[FIELD_WORDS]; // Starting with the codeword as the initial hash
    uint64_t combined[CONCAT_WORDS + CONCAT_WORDS]; // To hold codeword + hash or hash-only concatenations
    uint8_t digest[32];              // SHA3-256 output

    size_t path_offset = 0;          // Track current position in auth_path
    int current_idx = leaf_idx;      // Start at leaf index

    // Initialize current_hash with the leaf codeword
    memcpy(current_hash, leaf_codeword, FIELD_WORDS * sizeof(uint64_t));

    for (int layer = 0; layer <= 16; layer++) {
        if (layer == 0) {
            // Layer 0: Hash the leaf codeword with the sibling codeword element
            uint64_t *sibling_codeword = &auth_path[path_offset];
            path_offset += FIELD_WORDS;

            if (current_idx % 2 == 0) {
                memcpy(combined, current_hash, FIELD_WORDS * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS, sibling_codeword, FIELD_WORDS * sizeof(uint64_t));
            } else {
                memcpy(combined, sibling_codeword, FIELD_WORDS * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS, current_hash, FIELD_WORDS * sizeof(uint64_t));
            }

            // Hash the combined elements
            SHA3_host(digest, (uint8_t *)combined, 2 * FIELD_WORDS * sizeof(uint64_t), 256);
            memcpy(current_hash, digest, HASH_WORDS * sizeof(uint64_t));

        } else if (layer <= 11) {
            // Layers 1 to 11: Concatenate current and sibling's codeword and hash elements, then hash
            uint64_t *sibling_codeword = &auth_path[path_offset];
            path_offset += FIELD_WORDS;
            uint64_t *sibling_hash = &auth_path[path_offset];
            path_offset += HASH_WORDS;

            if (current_idx % 2 == 0) {
                memcpy(combined, current_hash, FIELD_WORDS * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS, sibling_codeword, FIELD_WORDS * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS + FIELD_WORDS, sibling_hash, HASH_WORDS * sizeof(uint64_t));
            } else {
                memcpy(combined, sibling_codeword, FIELD_WORDS * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS, sibling_hash, HASH_WORDS * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS + HASH_WORDS, current_hash, FIELD_WORDS * sizeof(uint64_t));
            }

            // Hash the combined elements for the next layer's current_hash
            SHA3_host(digest, (uint8_t *)combined, (FIELD_WORDS + FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t), 256);
            memcpy(current_hash, digest, HASH_WORDS * sizeof(uint64_t));

        } else if (layer <= 15) {
            // Layers 12 to 15: Only sibling hashes
            uint64_t *sibling_hash = &auth_path[path_offset];
            path_offset += HASH_WORDS;

            if (current_idx % 2 == 0) {
                memcpy(combined, current_hash, HASH_WORDS * sizeof(uint64_t));
                memcpy(combined + HASH_WORDS, sibling_hash, HASH_WORDS * sizeof(uint64_t));
            } else {
                memcpy(combined, sibling_hash, HASH_WORDS * sizeof(uint64_t));
                memcpy(combined + HASH_WORDS, current_hash, HASH_WORDS * sizeof(uint64_t));
            }

            // Hash the combined sibling hashes
            SHA3_host(digest, (uint8_t *)combined, 2 * HASH_WORDS * sizeof(uint64_t), 256);
            memcpy(current_hash, digest, HASH_WORDS * sizeof(uint64_t));
        }

        // Move up to the next layer’s parent index
        current_idx /= 2;
    }

    // After processing all layers, current_hash should contain the computed root
    return (memcmp(current_hash, expected_root, HASH_WORDS * sizeof(uint64_t)) == 0);
}
// int merkle_verify(
//     uint64_t *root,            // Expected Merkle root to verify against
//     size_t leaf_idx,           // Index of the leaf being verified
//     uint64_t **auth_path,      // Authentication path (sibling hashes) for the leaf
//     size_t proof_len,          // Length of the proof path
//     uint64_t *leaf             // Initial leaf (codeword element) to start the verification
// ) {
//     uint64_t computed_hash[HASH_WORDS];
//     memcpy(computed_hash, leaf, FIELD_WORDS * sizeof(uint64_t));

//     size_t current_idx = leaf_idx;

//     for (size_t i = 0; i < proof_len; ++i) {
//         uint64_t combined[CONCAT_WORDS];

//         if (i < proof_len - 5) {  // Middle layers with CONCAT_WORDS
//             if (current_idx % 2 == 0) {
//                 memcpy(combined, computed_hash, HASH_WORDS * sizeof(uint64_t));
//                 memcpy(combined + HASH_WORDS, auth_path[i], HASH_WORDS * sizeof(uint64_t));
//             } else {
//                 memcpy(combined, auth_path[i], HASH_WORDS * sizeof(uint64_t));
//                 memcpy(combined + HASH_WORDS, computed_hash, HASH_WORDS * sizeof(uint64_t));
//             }
//         } else {  // Last layers with HASH_WORDS
//             if (current_idx % 2 == 0) {
//                 memcpy(combined, computed_hash, HASH_WORDS * sizeof(uint64_t));
//                 memcpy(combined + HASH_WORDS, auth_path[i], HASH_WORDS * sizeof(uint64_t));
//             } else {
//                 memcpy(combined, auth_path[i], HASH_WORDS * sizeof(uint64_t));
//                 memcpy(combined + HASH_WORDS, computed_hash, HASH_WORDS * sizeof(uint64_t));
//             }
//         }

//         // Compute the hash for this layer
//         SHA3_host((uint8_t *)computed_hash, (uint8_t *)combined, CONCAT_WORDS * sizeof(uint64_t), 256);

//         // Move up to the parent node
//         current_idx /= 2;
//     }

//     // Final comparison with the provided root
//     return memcmp(computed_hash, root, HASH_WORDS * sizeof(uint64_t)) == 0;
// }