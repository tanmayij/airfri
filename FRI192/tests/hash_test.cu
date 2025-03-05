#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "../include/hash-host.cuh"  

#define FIELD_WORDS 4
#define HASH_WORDS 4
#define MAX_LAYERS 17  // Adjust based on your Merkle tree depth
#define MAX_NODES 131072  // Adjust based on your tree size

// Function prototype for SHA3 CPU computation
void SHA3_host(uint8_t *output, const uint8_t *input, size_t input_size, size_t bitSize);

// Structure to hold Merkle tree data
typedef struct {
    uint64_t ***tree;  // tree[layer][node][data]
    int layer_sizes[MAX_LAYERS];  // Stores number of nodes at each layer
} MerkleTree;

// Function to allocate and initialize Merkle tree structure
MerkleTree *initialize_merkle_tree() {
    MerkleTree *mtree = (MerkleTree *)malloc(sizeof(MerkleTree));
    if (!mtree) {
        fprintf(stderr, "Memory allocation failed for MerkleTree\n");
        return NULL;
    }

    mtree->tree = (uint64_t ***)malloc(MAX_LAYERS * sizeof(uint64_t **));
    if (!mtree->tree) {
        fprintf(stderr, "Memory allocation failed for tree structure\n");
        free(mtree);
        return NULL;
    }

    for (int i = 0; i < MAX_LAYERS; i++) {
        mtree->tree[i] = NULL;
        mtree->layer_sizes[i] = 0;
    }

    return mtree;
}

// Function to read and parse the Merkle tree from the given file format
MerkleTree *load_merkle_tree(const char *filename) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return NULL;
    }

    MerkleTree *mtree = initialize_merkle_tree();
    if (!mtree) {
        fclose(file);
        return NULL;
    }

    char line[512];
    int current_layer = -1, current_node = 0;
    
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "Layer")) {
            sscanf(line, "=== Layer %d", &current_layer);
            current_node = 0;

            if (mtree->tree[current_layer] == NULL) {
                mtree->tree[current_layer] = (uint64_t **)malloc(MAX_NODES * sizeof(uint64_t *));
            }
        } else if (strstr(line, "Node")) {
            int node_idx;
            uint64_t values[FIELD_WORDS + HASH_WORDS];

            int data_size = (current_layer == 0) ? FIELD_WORDS : (current_layer < 13) ? (FIELD_WORDS + HASH_WORDS) : HASH_WORDS;

            if (sscanf(line, "  Node %d: %lx %lx %lx %lx %lx %lx %lx %lx",
                       &node_idx, &values[0], &values[1], &values[2], &values[3],
                       &values[4], &values[5], &values[6], &values[7]) >= data_size) {

                mtree->tree[current_layer][node_idx] = (uint64_t *)malloc(data_size * sizeof(uint64_t));
                memcpy(mtree->tree[current_layer][node_idx], values, data_size * sizeof(uint64_t));

                mtree->layer_sizes[current_layer] = node_idx + 1;
            }
        }
    }

    fclose(file);
    return mtree;
}

// Function to verify SHA3 hashes in the Merkle tree
void verify_merkle_tree(MerkleTree *mtree, const char *output_filename) {
    FILE *output_file = fopen(output_filename, "w");
    if (!output_file) {
        fprintf(stderr, "Error opening output file: %s\n", output_filename);
        return;
    }

    for (int layer = 0; layer < MAX_LAYERS - 1; layer++) {
        int current_layer_size = mtree->layer_sizes[layer];
        int next_layer_size = mtree->layer_sizes[layer + 1];

        fprintf(output_file, "Layer %d Hash Verification:\n", layer);

        for (int i = 0; i < current_layer_size; i += 2) {
            if (i + 1 >= current_layer_size) break; // Avoid out-of-bounds

            uint64_t *combined;
            size_t combined_size;
            uint64_t computed_hash[HASH_WORDS];

            if (layer == 0) {
                // Layer 0: Codeword only
                combined_size = 2 * FIELD_WORDS;
                combined = (uint64_t *)malloc(combined_size * sizeof(uint64_t));
                memcpy(combined, mtree->tree[layer][i], FIELD_WORDS * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS, mtree->tree[layer][i + 1], FIELD_WORDS * sizeof(uint64_t));
            } else if (layer < 13) {
                // Layers 1-12: Codeword + Hash
                combined_size = 2 * (FIELD_WORDS + HASH_WORDS);
                combined = (uint64_t *)malloc(combined_size * sizeof(uint64_t));
                memcpy(combined, mtree->tree[layer][i], (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
                memcpy(combined + FIELD_WORDS + HASH_WORDS, mtree->tree[layer][i + 1], (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
            } else {
                // Layers 13-16: Hash only
                combined_size = 2 * HASH_WORDS;
                combined = (uint64_t *)malloc(combined_size * sizeof(uint64_t));
                memcpy(combined, mtree->tree[layer][i], HASH_WORDS * sizeof(uint64_t));
                memcpy(combined + HASH_WORDS, mtree->tree[layer][i + 1], HASH_WORDS * sizeof(uint64_t));
            }

            // Compute SHA3-256 hash
            SHA3_host((uint8_t *)computed_hash, (uint8_t *)combined, combined_size * sizeof(uint64_t), 256);
            free(combined);

            // Compare with the expected hash in the next layer
            int mismatch = memcmp(computed_hash, mtree->tree[layer + 1][i / 2], HASH_WORDS * sizeof(uint64_t));

            if (mismatch) {
                fprintf(output_file, "Mismatch at Layer %d Node %d\n", layer + 1, i / 2);
            }
        }
    }

    fclose(output_file);
}

// Free the allocated memory for the Merkle tree
void free_merkle_tree(MerkleTree *mtree) {
    for (int i = 0; i < MAX_LAYERS; i++) {
        if (mtree->tree[i]) {
            for (int j = 0; j < mtree->layer_sizes[i]; j++) {
                free(mtree->tree[i][j]);
            }
            free(mtree->tree[i]);
        }
    }
    free(mtree->tree);
    free(mtree);
}

// Main function
int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Usage: %s <input_file> <output_file>\n", argv[0]);
        return 1;
    }

    MerkleTree *mtree = load_merkle_tree(argv[1]);
    if (!mtree) {
        return 1;
    }

    verify_merkle_tree(mtree, argv[2]);

    free_merkle_tree(mtree);

    printf("Merkle tree verification completed. Mismatches (if any) are logged in %s\n", argv[2]);
    return 0;
}