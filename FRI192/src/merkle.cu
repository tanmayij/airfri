
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
const size_t CONCAT_WORDS = 8;  //FIELD_WORDS + HASH_WORDS

void print_field_merkle(const char *label, const uint64_t *field, int field_words) {
    printf("%s: ", label);
    for (int i = 0; i < field_words; i++) {
        printf("%016llx ", field[i]);
    }
    printf("\n");
}

//Helper function to hash data using SHA3-256
void hash_sha3_256(const uint64_t *data, size_t len, uint64_t *out) {
    EVP_MD_CTX *mdctx = EVP_MD_CTX_new();
    const EVP_MD *md = EVP_sha3_256();

    EVP_DigestInit_ex(mdctx, md, NULL);
    EVP_DigestUpdate(mdctx, data, len);
    EVP_DigestFinal_ex(mdctx, (unsigned char *)out, NULL);

    EVP_MD_CTX_free(mdctx);
}

//Function to print the byte array
void print_bytes(const uint64_t *bytes, size_t size) {
    for (size_t i = 0; i < size; i++) {
        printf("%016llx", bytes[i]);
        if (i < size - 1) {
            printf(" ");
        }
    }
    printf("\n");
}

//Function to compute the Merkle root from leaf hashes
void merkle_commit(uint64_t **leaf_hashes, size_t num_leaves, uint64_t *out) {
    //Assert to ensure memory is allocated correctly and no corruption occurs
    for (size_t i = 0; i < num_leaves; i++) {
        assert(leaf_hashes[i] != NULL && "leaf_hashes[i] is not properly allocated");
        //printf("leaf_hashes[%zu]: %p\n", i, (void*)leaf_hashes[i]);  //Debug to print pointers of leaf_hashes
    }

    assert((num_leaves & (num_leaves - 1)) == 0); //num_leaves must be a power of two

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

    //Compute the hash of the concatenated hashes
    hash_sha3_256(combined, 2 * HASH_SIZE, out);
    //printf("Computed root for current level\n");
}

//Function to check if a point (hash or codeword) has already been sent
bool is_sent(uint64_t *point, uint64_t **sent_points, size_t sent_count) {
    for (size_t i = 0; i < sent_count; i++) {
        if (memcmp(point, sent_points[i], HASH_SIZE) == 0) {
            return true;
        }
    }
    return false;
}
//can have merkle_open not worry about the tree data structure, but compute the auth path 
//to do that, we need a data structure that handles the hashes only? then we could pick an index, get the corresponding codeword element and the hash element for that index.

void merkle_open(uint64_t **auth_path, int leaf_idx, size_t *proof_len, uint64_t ***tree, int layer) {
    int current_index = leaf_idx;
    *proof_len = 0;

    for (int i = layer; i < 17; i++) {
        int sibling_index = (current_index % 2 == 0) ? current_index + 1 : current_index - 1;
        size_t sibling_size;

        if (i == 0) {
            //For layer 0, only copy the sibling (FIELD_WORDS size)
            sibling_size = FIELD_WORDS;
            auth_path[*proof_len] = (uint64_t *)malloc(sibling_size * sizeof(uint64_t));
            if (!auth_path[*proof_len]) {
                fprintf(stderr, "Memory allocation failed for auth_path at proof_len %zu\n", *proof_len);
                exit(1);
            }
            memcpy(auth_path[*proof_len], tree[i][sibling_index], sibling_size * sizeof(uint64_t));

        } else if (i > 0 && i < 13) {
            //For layers 1 to 12, include the codeword element (FIELD_WORDS) and sibling
            sibling_size = FIELD_WORDS + HASH_WORDS;
            auth_path[*proof_len] = (uint64_t *)malloc((FIELD_WORDS + sibling_size) * sizeof(uint64_t));
            if (!auth_path[*proof_len]) {
                fprintf(stderr, "Memory allocation failed for auth_path at proof_len %zu\n", *proof_len);
                exit(1);
            }
            //Copy the codeword element from the current index (skipping the hash part)

            memcpy(auth_path[*proof_len], tree[i][current_index], FIELD_WORDS * sizeof(uint64_t));
            memcpy(auth_path[*proof_len] + FIELD_WORDS, tree[i][sibling_index], sibling_size * sizeof(uint64_t));

        } else {
            //For layers 13 to 16, only include the sibling (HASH_WORDS size)
            sibling_size = HASH_WORDS;
            auth_path[*proof_len] = (uint64_t *)malloc(sibling_size * sizeof(uint64_t));
            if (!auth_path[*proof_len]) {
                fprintf(stderr, "Memory allocation failed for auth_path at proof_len %zu\n", *proof_len);
                exit(1);
            }
            memcpy(auth_path[*proof_len], tree[i][sibling_index], sibling_size * sizeof(uint64_t));
        }

        (*proof_len)++;
        current_index /= 2;  //Move up the tree for the next layer
        //printf("current index is: %d\n", current_index);
    }
}

//merkle_verify should do the following:
//step 1: when i = 0, combine the leaf and the first element (which will be codeword) in the auth path.
//step 2: when i > 0 and i < 13, combine the codeword and the hash in the auth path. when leaf_idx is odd, add computed hash to the end of the combined array. when leaf_idx is even, add computed hash starting the [FIELD_WORDS + 1] index of the combined array.
//step 3: when i >= 13, combine the two hashes in the auth path. same as step 2, add computed hash to the end of the combined array if leaf_idx is odd, otherwise add computed hash starting the [HASH_WORDS] index of the combined array.
//step 4: compute the new hash using the combined array.
//step 5: compare the computed hash with the root hash. return 1 if they are equal, 0 otherwise.

//Function to verify a Merkle proof
// int merkle_verify(
//     uint64_t *root,           
//     size_t leaf_idx,          
//     uint64_t **auth_path,     
//     size_t proof_len,         
//     uint64_t *leaf,           
//     int layer
// ) {
//     // Allocate memory for current_hash
//     uint64_t *current_hash = (uint64_t *)malloc(HASH_WORDS * sizeof(uint64_t));
//     if (current_hash == NULL) {
//         fprintf(stderr, "Memory allocation failed for current_hash\n");
//         exit(1);
//     }

//     // Allocate memory for new_hash
//     uint64_t *new_hash = (uint64_t *)calloc(HASH_WORDS, sizeof(uint64_t));
//     if (new_hash == NULL) {
//         fprintf(stderr, "Memory allocation failed for new_hash\n");
//         free(current_hash);
//         exit(1);
//     }

//     memcpy(current_hash, leaf, FIELD_WORDS * sizeof(uint64_t));

//     for (size_t i = layer; i < proof_len; i++) {
//         uint64_t *combined;
//         size_t combined_size;

//         int is_even = (leaf_idx % 2 == 0);

//         if (i == 0) {
//             // **Layer 0**: Only combine 2 * FIELD_WORDS
//             combined_size = 2 * FIELD_WORDS;
//             combined = (uint64_t *)malloc(combined_size * sizeof(uint64_t));

//             if (combined == NULL) {
//                 fprintf(stderr, "Memory allocation failed for combined\n");
//                 free(current_hash);
//                 free(new_hash);
//                 exit(1);
//             }

//             memcpy(combined, current_hash, FIELD_WORDS * sizeof(uint64_t));
//             memcpy(combined + FIELD_WORDS, auth_path[i], FIELD_WORDS * sizeof(uint64_t));

//         } else if (i < 13) {
//             // **Layers 1-12**: Lone codeword + hash handling
//             if (is_even) {
//                 // Even index: lone codeword + computed hash + (codeword + hashword)
//                 combined_size = 2 * FIELD_WORDS + HASH_WORDS;
//                 combined = (uint64_t *)malloc(combined_size * sizeof(uint64_t));

//                 if (combined == NULL) {
//                     fprintf(stderr, "Memory allocation failed for combined\n");
//                     free(current_hash);
//                     free(new_hash);
//                     exit(1);
//                 }

//                 memcpy(combined, auth_path[i], FIELD_WORDS * sizeof(uint64_t));
//                 memcpy(combined + FIELD_WORDS, current_hash, HASH_WORDS * sizeof(uint64_t));
//                 memcpy(combined + FIELD_WORDS + HASH_WORDS, auth_path[i] + FIELD_WORDS, (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
//             } 
//             else {
//                 // Odd index: (codeword + hashword) + lone codeword + computed hash
//                 combined_size = 2 * FIELD_WORDS + HASH_WORDS;
//                 combined = (uint64_t *)malloc(combined_size * sizeof(uint64_t));

//                 if (combined == NULL) {
//                     fprintf(stderr, "Memory allocation failed for combined\n");
//                     free(current_hash);
//                     free(new_hash);
//                     exit(1);
//                 }

//                 memcpy(combined, auth_path[i] + FIELD_WORDS, (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
//                 memcpy(combined + FIELD_WORDS + HASH_WORDS, auth_path[i], FIELD_WORDS * sizeof(uint64_t));
//                 memcpy(combined + FIELD_WORDS, current_hash, HASH_WORDS * sizeof(uint64_t));
//             }
//         } 
//         else {
//             // **Layers 13-16**: Standard `HASH_WORDS + HASH_WORDS`
//             combined_size = 2 * HASH_WORDS;
//             combined = (uint64_t *)malloc(combined_size * sizeof(uint64_t));

//             if (combined == NULL) {
//                 fprintf(stderr, "Memory allocation failed for combined\n");
//                 free(current_hash);
//                 free(new_hash);
//                 exit(1);
//             }

//             if (is_even) {
//                 memcpy(combined, current_hash, HASH_WORDS * sizeof(uint64_t));
//                 memcpy(combined + HASH_WORDS, auth_path[i], HASH_WORDS * sizeof(uint64_t));
//             } else {
//                 memcpy(combined, auth_path[i], HASH_WORDS * sizeof(uint64_t));
//                 memcpy(combined + HASH_WORDS, current_hash, HASH_WORDS * sizeof(uint64_t));
//             }
//         }

//         // Compute new hash
//         SHA3_host((uint8_t *)new_hash, (uint8_t *)combined, combined_size * sizeof(uint64_t), 256);
        
//         memcpy(current_hash, new_hash, HASH_WORDS * sizeof(uint64_t));
//         memset(new_hash, 0, HASH_WORDS * sizeof(uint64_t));

//         // Debugging print
//         printf("Layer %zu Debug:\n", i);
//         printf("Combined Input: ");
//         for (size_t j = 0; j < combined_size; j++) {
//             printf("%016lx ", combined[j]);
//         }
//         printf("\nNew Current Hash: ");
//         for (int j = 0; j < HASH_WORDS; j++) {
//             printf("%016lx ", current_hash[j]);
//         }
//         printf("\n");

//         free(combined);
//         leaf_idx /= 2;  // Move up the tree for the next layer
//     }

//     printf("Computed Merkle Root: ");
//     for (int i = 0; i < HASH_WORDS; i++) {
//         printf("%016lx ", current_hash[i]);
//     }
//     printf("\n");

//     int result = memcmp(root, current_hash, HASH_WORDS * sizeof(uint64_t)) == 0;

//     free(current_hash);
//     free(new_hash);

//     return result;
// }

int merkle_verify(
    uint64_t *root,
    size_t leaf_idx,
    uint64_t **auth_path,
    size_t proof_len,
    uint64_t *leaf,
    int layer,
    int protocol_layer // New parameter to handle FRI round
) {
    uint64_t *current_hash = (uint64_t *)malloc(HASH_WORDS * sizeof(uint64_t));
    if (current_hash == NULL) {
        fprintf(stderr, "Memory allocation failed for current_hash\n");
        exit(1);
    }

    uint64_t *new_hash = (uint64_t *)calloc(HASH_WORDS, sizeof(uint64_t));
    if (new_hash == NULL) {
        fprintf(stderr, "Memory allocation failed for new_hash\n");
        free(current_hash);
        exit(1);
    }

    if (protocol_layer == 0) {
        // First round: use the given leaf value
        memcpy(current_hash, leaf, FIELD_WORDS * sizeof(uint64_t));
    } else {
        // Second and subsequent rounds: use the auth_path[layer] value
        memcpy(current_hash, auth_path[layer], HASH_WORDS * sizeof(uint64_t));
        layer++; // Skip the first auth_path since it's already used
    }

    for (size_t i = layer; i < proof_len; i++) {
        size_t combined_size;
        uint64_t *combined;

        int is_even = (leaf_idx % 2 == 0);

        if (i == 0 && protocol_layer == 0) {
            combined_size = 2 * FIELD_WORDS;
            combined = (uint64_t *)malloc(combined_size * sizeof(uint64_t));

            if (!combined) {
                fprintf(stderr, "Memory allocation failed\n");
                free(current_hash);
                free(new_hash);
                exit(1);
            }

            if (is_even) {
                memcpy(combined, current_hash, FIELD_WORDS * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS, auth_path[i], FIELD_WORDS * sizeof(uint64_t));
            } else {
                memcpy(combined, auth_path[i], FIELD_WORDS * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS, current_hash, FIELD_WORDS * sizeof(uint64_t));
            }
        }

        else if (i < 13) {
            combined_size = 2 * (FIELD_WORDS + HASH_WORDS);
            combined = (uint64_t *)malloc(combined_size * sizeof(uint64_t));

            if (!combined) {
                fprintf(stderr, "Memory allocation failed\n");
                free(current_hash);
                free(new_hash);
                exit(1);
            }

            if (is_even) {
                // Even index: Lone Codeword + Computed Hash + (Codeword + Hashword)
                memcpy(combined, auth_path[i], FIELD_WORDS * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS, current_hash, HASH_WORDS * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS + HASH_WORDS, auth_path[i] + FIELD_WORDS, (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
            } else {
                // Odd index: (Codeword + Hashword) + Lone Codeword + Computed Hash
                memcpy(combined, auth_path[i] + FIELD_WORDS, (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS + HASH_WORDS, auth_path[i], FIELD_WORDS * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS + HASH_WORDS + FIELD_WORDS, current_hash, HASH_WORDS * sizeof(uint64_t));
            }
        }

        else {
            combined_size = 2 * HASH_WORDS;
            combined = (uint64_t *)malloc(combined_size * sizeof(uint64_t));

            if (!combined) {
                fprintf(stderr, "Memory allocation failed\n");
                free(current_hash);
                free(new_hash);
                exit(1);
            }

            if (is_even) {
                memcpy(combined, current_hash, HASH_WORDS * sizeof(uint64_t));
                memcpy(combined + HASH_WORDS, auth_path[i], HASH_WORDS * sizeof(uint64_t));
            } else {
                memcpy(combined, auth_path[i], HASH_WORDS * sizeof(uint64_t));
                memcpy(combined + HASH_WORDS, current_hash, HASH_WORDS * sizeof(uint64_t));
            }
        }

        // Hash the combined data
        SHA3_host((uint8_t *)new_hash, (uint8_t *)combined, combined_size * sizeof(uint64_t), 256);
        memcpy(current_hash, new_hash, HASH_WORDS * sizeof(uint64_t));
        memset(new_hash, 0, HASH_WORDS * sizeof(uint64_t));

        printf("Layer %zu Debug:\n", i);
        printf("Combined Input: ");
        for (size_t j = 0; j < combined_size; j++) {
            printf("%016lx ", combined[j]);
        }
        printf("\nNew Current Hash: ");
        for (int j = 0; j < HASH_WORDS; j++) {
            printf("%016lx ", current_hash[j]);
        }
        printf("\n");

        free(combined);
        leaf_idx /= 2;
    }

    printf("Computed Merkle Root: ");
    for (int i = 0; i < HASH_WORDS; i++) {
        printf("%016lx ", current_hash[i]);
    }
    printf("\n");

    int result = memcmp(root, current_hash, HASH_WORDS * sizeof(uint64_t)) == 0;

    free(current_hash);
    free(new_hash);

    return result;
}