// #include <openssl/evp.h>
// #include <stdint.h>
// #include <stdlib.h>
// #include <stdbool.h>
// #include <string.h>
// #include <assert.h>
// #include "../include/merkle.cuh"
// #include "../include/hash-host.cuh"
// #include <stdio.h>

// void print_field_merkle(const char *label, const uint64_t *field, int field_words) {
//     printf("%s: ", label);
//     for (int i = 0; i < field_words; i++) {
//         printf("%016llx ", field[i]);
//     }
//     printf("\n");
// }

// // Helper function to hash data using SHA3-256
// void hash_sha3_256(const uint64_t *data, size_t len, uint64_t *out) {
//     EVP_MD_CTX *mdctx = EVP_MD_CTX_new();
//     const EVP_MD *md = EVP_sha3_256();

//     EVP_DigestInit_ex(mdctx, md, NULL);
//     EVP_DigestUpdate(mdctx, data, len);
//     EVP_DigestFinal_ex(mdctx, (unsigned char *)out, NULL);

//     EVP_MD_CTX_free(mdctx);
// }

// // Function to print the byte array
// void print_bytes(const uint64_t *bytes, size_t size) {
//     for (size_t i = 0; i < size; i++) {
//         printf("%016llx", bytes[i]);
//         if (i < size - 1) {
//             printf(" ");
//         }
//     }
//     printf("\n");
// }

// // Function to compute the Merkle root from leaf hashes
// void merkle_commit(uint64_t **leaf_hashes, size_t num_leaves, uint64_t *out) {
//     // Assert to ensure memory is allocated correctly and no corruption occurs
//     for (size_t i = 0; i < num_leaves; i++) {
//         assert(leaf_hashes[i] != NULL && "leaf_hashes[i] is not properly allocated");
//         //printf("leaf_hashes[%zu]: %p\n", i, (void*)leaf_hashes[i]);  // Debug to print pointers of leaf_hashes
//     }

//     assert((num_leaves & (num_leaves - 1)) == 0); // num_leaves must be a power of two

//     if (num_leaves == 1) {
//         //printf("Computing Merkle root for single leaf\n");
//         memcpy(out, leaf_hashes[0], HASH_SIZE);
//         return;
//     }

//     size_t half = num_leaves / 2;
//     uint64_t left_root[HASH_WORDS];
//     uint64_t right_root[HASH_WORDS];

//     //printf("Recursing left on %zu leaves\n", half);
//     merkle_commit(leaf_hashes, half, left_root);
//     //printf("Recursing right on %zu leaves\n", half);
//     merkle_commit(leaf_hashes + half, half, right_root);

//     uint64_t combined[2 * HASH_WORDS];
//     memcpy(combined, left_root, HASH_SIZE);
//     memcpy(combined + HASH_WORDS, right_root, HASH_SIZE);

//     // Compute the hash of the concatenated hashes
//     hash_sha3_256(combined, 2 * HASH_SIZE, out);
//     //printf("Computed root for current level\n");
// }

// // Function to check if a point (hash or codeword) has already been sent
// bool is_sent(uint64_t *point, uint64_t **sent_points, size_t sent_count) {
//     for (size_t i = 0; i < sent_count; i++) {
//         if (memcmp(point, sent_points[i], HASH_SIZE) == 0) {
//             return true;
//         }
//     }
//     return false;
// }

// void merkle_open(
//     uint64_t ***codewords,      // 3D array for codewords at each layer
//     size_t num_leaves,          // Number of leaves in the tree (length of initial codewords)
//     size_t leaf_idx,            // Index of the leaf element
//     uint64_t **proof_path,      // Output array for proof path
//     size_t *proof_len,          // Output proof path length
//     size_t field_words          // Number of words in each field element
// ) {
//     *proof_len = 0;
//     size_t index = leaf_idx;
//     size_t codeword_layers = 12; //always 12 layers of codewords 
//     uint64_t combined[2 * field_words];  // Buffer to combine elements

//     // Traverse up the tree through codeword layers
//     for (size_t layer = 0; layer < codeword_layers; layer++) {
//         size_t sibling_index = (index % 2 == 0) ? index + 1 : index - 1;

//         // Store the current codeword (f_xx) in the proof path
//         proof_path[*proof_len] = (uint64_t *)malloc(field_words * sizeof(uint64_t));
//         memcpy(proof_path[*proof_len], codewords[layer][index], field_words * sizeof(uint64_t));
//         (*proof_len)++;

//         // Combine current codeword with sibling codeword for hash storage
//         memcpy(combined, codewords[layer][index], field_words * sizeof(uint64_t));
//         memcpy(combined + field_words, codewords[layer][sibling_index], field_words * sizeof(uint64_t));

//         // Hash the combined codeword elements and store in proof path
//         proof_path[*proof_len] = (uint64_t *)malloc(HASH_WORDS * sizeof(uint64_t));
//         SHA3_host((uint8_t *)proof_path[*proof_len], (uint8_t *)combined, 2 * field_words * sizeof(uint64_t), 256);
//         (*proof_len)++;

//         index /= 2;  // Move to the parent node
//     }
//     int num_layers = (int)log2(num_leaves);
//     // Continue up the tree with hashes only beyond codeword layers
//     for (size_t layer = codeword_layers; layer < num_layers; layer++) {  
//         size_t sibling_index = (index % 2 == 0) ? index + 1 : index - 1;
//         uint64_t combined_hash[2 * HASH_WORDS];

//         // Combine current hash with sibling hash
//         if (index % 2 == 0) {
//             memcpy(combined_hash, proof_path[*proof_len - 1], HASH_WORDS * sizeof(uint64_t));
//             memcpy(combined_hash + HASH_WORDS, codewords[codeword_layers - 1][sibling_index], HASH_WORDS * sizeof(uint64_t));
//         } else {
//             memcpy(combined_hash, codewords[codeword_layers - 1][sibling_index], HASH_WORDS * sizeof(uint64_t));
//             memcpy(combined_hash + HASH_WORDS, proof_path[*proof_len - 1], HASH_WORDS * sizeof(uint64_t));
//         }

//         // Hash the combined array and add it to the proof path
//         proof_path[*proof_len] = (uint64_t *)malloc(HASH_WORDS * sizeof(uint64_t));
//         SHA3_host((uint8_t *)proof_path[*proof_len], (uint8_t *)combined, 2 * HASH_WORDS * sizeof(uint64_t), 256);
//         (*proof_len)++;

//         index /= 2;  // Move up the tree
//     }
// }

// bool merkle_verify(
//     uint64_t *root_hash,          // Known root hash of the Merkle tree
//     uint64_t *leaf,               // Initial codeword element (leaf node)
//     size_t leaf_idx,              // Index of the leaf element
//     uint64_t **proof_path,        // Proof path to verify
//     size_t proof_len,             // Length of the proof path
//     size_t field_words            // Number of words in each field element
// ) {
//     size_t index = leaf_idx;
//     uint64_t *current_hash = (uint64_t *)malloc(HASH_WORDS * sizeof(uint64_t));

//     // Hash the initial codeword (leaf) if field_words > 1
//     if (field_words > 1) {
//         SHA3_host((uint8_t *)current_hash, (uint8_t *)leaf, field_words * sizeof(uint64_t), 256);
//     } else {
//         // If field_words == 1, directly use the leaf value
//         memcpy(current_hash, leaf, field_words * sizeof(uint64_t));
//     }

//     // Traverse the proof path
//     for (size_t i = 0; i < proof_len; i++) {
//         uint64_t combined[2 * field_words];
        
//         // For layers that combine codewords and hashes
//         if (i % 2 == 0 && i < proof_len - 1) {
//             // Combine codeword and hash
//             if (index % 2 == 0) {
//                 memcpy(combined, current_hash, field_words * sizeof(uint64_t));
//                 memcpy(combined + field_words, proof_path[i + 1], field_words * sizeof(uint64_t));
//             } else {
//                 memcpy(combined, proof_path[i + 1], field_words * sizeof(uint64_t));
//                 memcpy(combined + field_words, current_hash, field_words * sizeof(uint64_t));
//             }
//             // Hash this combination and update current_hash
//             SHA3_host((uint8_t *)current_hash, (uint8_t *)combined, 2 * field_words * sizeof(uint64_t), 256);
//             i++;  // Skip next proof_path element since we've already used it
//         } else {
//             // When only hashes are left, just combine them as usual
//             uint64_t combined_hash[2 * HASH_WORDS];
//             if (index % 2 == 0) {
//                 memcpy(combined_hash, current_hash, HASH_WORDS * sizeof(uint64_t));
//                 memcpy(combined_hash + HASH_WORDS, proof_path[i], HASH_WORDS * sizeof(uint64_t));
//             } else {
//                 memcpy(combined_hash, proof_path[i], HASH_WORDS * sizeof(uint64_t));
//                 memcpy(combined_hash + HASH_WORDS, current_hash, HASH_WORDS * sizeof(uint64_t));
//             }
//             SHA3_host((uint8_t *)current_hash, (uint8_t *)combined_hash, 2 * HASH_WORDS * sizeof(uint64_t), 256);
//         }
        
//         index /= 2;  // Move up to the parent index

//         // Debug print to track current hash after each proof step
//         printf("Step %zu, Computed Hash: ", i);
//         for (size_t j = 0; j < HASH_WORDS; j++) {
//             printf("%016lx ", current_hash[j]);
//         }
//         printf("\n");
//     }

//     // Check if the computed hash matches the root hash
//     bool result = (memcmp(current_hash, root_hash, HASH_WORDS * sizeof(uint64_t)) == 0);

//     // Final debug output for comparison with the root hash
//     printf("Final Computed Hash: ");
//     for (size_t i = 0; i < HASH_WORDS; i++) {
//         printf("%016lx ", current_hash[i]);
//     }
//     printf("\n");

//     printf("Root Hash: ");
//     for (size_t i = 0; i < HASH_WORDS; i++) {
//         printf("%016lx ", root_hash[i]);
//     }
//     printf("\n");

//     free(current_hash);
//     return result;
// }
// // bool merkle_verify(
// //     uint64_t *root_hash,          // Known root hash of the Merkle tree
// //     uint64_t *leaf,               // Initial codeword element (leaf node)
// //     size_t leaf_idx,              // Index of the leaf element
// //     uint64_t **proof_path,        // Proof path to verify
// //     size_t proof_len,             // Length of the proof path
// //     size_t field_words            // Number of words in each field element
// // ) {
// //     size_t index = leaf_idx;
// //     uint64_t *current_hash = (uint64_t *)malloc(HASH_WORDS * sizeof(uint64_t));

// //     // Hash the initial codeword (leaf) if field_words > 1
// //     if (field_words > 1) {
// //         SHA3_host((uint8_t *)current_hash, (uint8_t *)leaf, field_words * sizeof(uint64_t), 256);
// //     } else {
// //         // If field_words == 1, directly use the leaf value
// //         memcpy(current_hash, leaf, field_words * sizeof(uint64_t));
// //     }

// //     // Traverse the proof path
// //     for (size_t i = 0; i < proof_len; i++) {
// //         uint64_t combined[2 * field_words];
        
// //         // Odd levels may include codewords in the path
// //         if (i < proof_len && i % 2 == 0 && i < proof_len - 1) {
// //             // Combine codeword and hash
// //             if (index % 2 == 0) {
// //                 memcpy(combined, current_hash, field_words * sizeof(uint64_t));
// //                 memcpy(combined + field_words, proof_path[i + 1], field_words * sizeof(uint64_t));
// //             } else {
// //                 memcpy(combined, proof_path[i + 1], field_words * sizeof(uint64_t));
// //                 memcpy(combined + field_words, current_hash, field_words * sizeof(uint64_t));
// //             }
// //             // Hash this combination and update current_hash
// //             SHA3_host((uint8_t *)current_hash, (uint8_t *)combined, 2 * field_words * sizeof(uint64_t), 256);
// //             i++;  // Skip next proof_path since we've already used it
// //         } else {
// //             // When only hashes are left, just combine them as usual
// //             uint64_t combined_hash[2 * HASH_WORDS];
// //             if (index % 2 == 0) {
// //                 memcpy(combined_hash, current_hash, HASH_WORDS * sizeof(uint64_t));
// //                 memcpy(combined_hash + HASH_WORDS, proof_path[i], HASH_WORDS * sizeof(uint64_t));
// //             } else {
// //                 memcpy(combined_hash, proof_path[i], HASH_WORDS * sizeof(uint64_t));
// //                 memcpy(combined_hash + HASH_WORDS, current_hash, HASH_WORDS * sizeof(uint64_t));
// //             }
// //             SHA3_host((uint8_t *)current_hash, (uint8_t *)combined_hash, 2 * HASH_WORDS * sizeof(uint64_t), 256);
// //         }
        
// //         index /= 2;  // Move up to the parent index
// //     }

// //     // Check if the computed hash matches the root hash
// //     bool result = (memcmp(current_hash, root_hash, HASH_WORDS * sizeof(uint64_t)) == 0);
// //     printf("Computed Hash: ");
// //     for (size_t i = 0; i < HASH_WORDS; i++) {
// //         printf("%016lx ", current_hash[i]);
// //     }
// //     printf("\n");
// //     printf("Actual Hash: ");
// //     for (size_t i = 0; i < HASH_WORDS; i++) {
// //         printf("%016lx ", root_hash[i]);
// //     }
// //     free(current_hash);
// //     return result;
// // }


// // // Example main function to test the Merkle functions
// // int main() {
// //     // Example data: an array of data elements (integers here for simplicity)
// //     int data[] = {2, 3, 4, 5};
// //     size_t num_leaves = sizeof(data) / sizeof(data[0]);

// //     // Convert integers to strings and then to hashes
// //     uint64_t *leaf_hashes[num_leaves];
// //     for (size_t i = 0; i < num_leaves; i++) {
// //         // Convert integer to string
// //         char str[12]; // Large enough to hold any 32-bit integer
// //         snprintf(str, sizeof(str), "%d", data[i]);

// //         // Allocate memory for the hash
// //         leaf_hashes[i] = (uint64_t *)malloc(HASH_SIZE);
// //         if (leaf_hashes[i] == NULL) {
// //             fprintf(stderr, "Memory allocation failed for leaf hash %zu\n", i);
// //             return 1;
// //         }

// //         // Compute the hash
// //         hash_sha3_256((uint64_t *)str, strlen(str), leaf_hashes[i]);
// //     }

// //     // Print each hash in the leaf_hashes array
// //     for (size_t i = 0; i < num_leaves; i++) {
// //         printf("Hash of leaf %zu: ", i);
// //         print_bytes(leaf_hashes[i], HASH_WORDS);
// //     }

// //     // Compute the Merkle root
// //     uint64_t merkle_root[HASH_WORDS];
// //     merkle_commit(leaf_hashes, num_leaves, merkle_root);

// //     // Print the Merkle root
// //     printf("Merkle Root: ");
// //     print_bytes(merkle_root, HASH_WORDS);

// //     // Generate a proof for a specific leaf (index 0 in this example)
// //     uint64_t *proof_path[num_leaves];
// //     size_t proof_len = 0;
// //     merkle_open(leaf_hashes, num_leaves, 0, proof_path, &proof_len);

// //     // Print the proof path
// //     printf("Proof Path: ");
// //     for (size_t i = 0; i < proof_len; i++) {
// //         print_bytes(proof_path[i], HASH_WORDS);
// //         if (i < proof_len - 1) printf(" ");
// //     }
// //     printf("\n");

// //     // Verify the proof
// //     int is_valid = merkle_verify(merkle_root, 0, proof_path, proof_len, leaf_hashes[0]);
// //     printf("Verification result: %s\n", is_valid ? "valid" : "invalid");

// //     // Clean up allocated memory
// //     for (size_t i = 0; i < num_leaves; i++) {
// //         free(leaf_hashes[i]);
// //     }

// //     // // Free proof path memory
// //     // for (size_t i = 0; i < proof_len; i++) {
// //     //     free(proof_path[i]);
// //     // }

// //     return 0;
// // }

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

void merkle_open(uint64_t **auth_path, int leaf_idx, size_t *proof_len, uint64_t *host_concatenated_tree) {
    int current_index = leaf_idx;
    *proof_len = 0;  // Initialize proof_len to 0 and increment as we add layers

    // Traverse up the tree and collect sibling hashes
    while (current_index > 0) {
        int sibling_index;

        // Determine sibling index: if current index is even, sibling is current_index + 1; if odd, it's current_index - 1
        if (current_index % 2 == 0) {
            sibling_index = current_index + 1;
        } else {
            sibling_index = current_index - 1;
        }

        // Allocate space for the sibling hash in auth_path
        auth_path[*proof_len] = (uint64_t *)malloc(CONCAT_WORDS * sizeof(uint64_t));

        // Copy the sibling hash from host_concatenated_tree to auth_path for the current layer
        memcpy(auth_path[*proof_len], 
               host_concatenated_tree + sibling_index * CONCAT_WORDS, 
               CONCAT_WORDS * sizeof(uint64_t));

        // Increment proof_len since we're adding one layer to the proof path
        (*proof_len)++;

        // Move up one level in the tree (integer division to find parent index)
        current_index /= 2;
    }
}
// void merkle_open(
//     uint64_t **auth_path,          // 2D array to store the sibling hashes at each layer
//     int leaf_idx,                  // Index of the leaf for which the path is generated
//     size_t *proof_len,             // Will store the length of the proof path
//     uint64_t *tree_initial_leaf,   // Array holding the initial codeword elements
//     uint64_t *tree_concat_words,   // Array for nodes that are CONCAT_WORDS-sized (codeword || hash)
//     uint64_t *tree_hash_words,     // Array for nodes that are only HASH_WORDS-sized
//     size_t num_layers              // Number of layers in the tree
// ) {
//     int current_index = leaf_idx;
//     *proof_len = 0;

//     // Traverse up the tree, layer by layer
//     for (size_t layer = 0; layer < num_layers; ++layer) {
//         int sibling_index = (current_index % 2 == 0) ? current_index + 1 : current_index - 1;

//         // Select the correct tree array based on the layer type
//         if (layer == 0) {  // First layer with codewords only
//             auth_path[*proof_len] = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
//             memcpy(auth_path[*proof_len], &tree_initial_leaf[sibling_index * FIELD_WORDS], FIELD_WORDS * sizeof(uint64_t));
//         } else if (layer < num_layers - 5) {  // Middle layers with CONCAT_WORDS (codeword || hash)
//             auth_path[*proof_len] = (uint64_t *)malloc(CONCAT_WORDS * sizeof(uint64_t));
//             memcpy(auth_path[*proof_len], &tree_concat_words[sibling_index * CONCAT_WORDS], CONCAT_WORDS * sizeof(uint64_t));
//         } else {  // Last layers with HASH_WORDS only
//             auth_path[*proof_len] = (uint64_t *)malloc(HASH_WORDS * sizeof(uint64_t));
//             memcpy(auth_path[*proof_len], &tree_hash_words[(sibling_index - (1 << (num_layers - 5))) * HASH_WORDS], HASH_WORDS * sizeof(uint64_t));
//         }

//         (*proof_len)++;
//         current_index /= 2;  // Move up to the parent node
//     }
// }

// int merkle_verify(uint64_t *root, size_t index, uint64_t **auth_path, size_t proof_len, uint64_t *leaf) {
//     uint64_t computed_hash[HASH_WORDS];
//     memcpy(computed_hash, leaf, FIELD_WORDS * sizeof(uint64_t));

//     for (size_t i = 0; i < proof_len; ++i) {
//         uint64_t combined[CONCAT_WORDS];

//         if (index % 2 == 0) {
//             // Left node: concatenate computed hash with the proof path element
//             memcpy(combined, computed_hash, HASH_WORDS * sizeof(uint64_t));
//             memcpy(combined + HASH_WORDS, auth_path[i], HASH_WORDS * sizeof(uint64_t));
//         } else {
//             // Right node: concatenate proof path element with the computed hash
//             memcpy(combined, auth_path[i], HASH_WORDS * sizeof(uint64_t));
//             memcpy(combined + HASH_WORDS, computed_hash, HASH_WORDS * sizeof(uint64_t));
//         }

//         // Recompute hash
//         SHA3_host((uint8_t *)computed_hash, (uint8_t *)combined, CONCAT_WORDS * sizeof(uint64_t), 256);
//         index /= 2;

//         // Print each computed hash for debugging
//         printf("Layer %zu: Computed hash = ", i);
//         for (size_t j = 0; j < HASH_WORDS; ++j) {
//             printf("%016lx ", computed_hash[j]);
//         }
//         printf("\n");
//     }

//     // Final comparison with the root
//     return memcmp(computed_hash, root, HASH_WORDS * sizeof(uint64_t)) == 0;
// }
int merkle_verify(
    uint64_t *root,            // Expected Merkle root to verify against
    size_t leaf_idx,           // Index of the leaf being verified
    uint64_t **auth_path,      // Authentication path (sibling hashes) for the leaf
    size_t proof_len,          // Length of the proof path
    uint64_t *leaf             // Initial leaf (codeword element) to start the verification
) {
    uint64_t computed_hash[HASH_WORDS];
    memcpy(computed_hash, leaf, FIELD_WORDS * sizeof(uint64_t));

    size_t current_idx = leaf_idx;

    for (size_t i = 0; i < proof_len; ++i) {
        uint64_t combined[CONCAT_WORDS];

        if (i < proof_len - 5) {  // Middle layers with CONCAT_WORDS
            if (current_idx % 2 == 0) {
                memcpy(combined, computed_hash, HASH_WORDS * sizeof(uint64_t));
                memcpy(combined + HASH_WORDS, auth_path[i], HASH_WORDS * sizeof(uint64_t));
            } else {
                memcpy(combined, auth_path[i], HASH_WORDS * sizeof(uint64_t));
                memcpy(combined + HASH_WORDS, computed_hash, HASH_WORDS * sizeof(uint64_t));
            }
        } else {  // Last layers with HASH_WORDS
            if (current_idx % 2 == 0) {
                memcpy(combined, computed_hash, HASH_WORDS * sizeof(uint64_t));
                memcpy(combined + HASH_WORDS, auth_path[i], HASH_WORDS * sizeof(uint64_t));
            } else {
                memcpy(combined, auth_path[i], HASH_WORDS * sizeof(uint64_t));
                memcpy(combined + HASH_WORDS, computed_hash, HASH_WORDS * sizeof(uint64_t));
            }
        }

        // Compute the hash for this layer
        SHA3_host((uint8_t *)computed_hash, (uint8_t *)combined, CONCAT_WORDS * sizeof(uint64_t), 256);

        // Move up to the parent node
        current_idx /= 2;
    }

    // Final comparison with the provided root
    return memcmp(computed_hash, root, HASH_WORDS * sizeof(uint64_t)) == 0;
}