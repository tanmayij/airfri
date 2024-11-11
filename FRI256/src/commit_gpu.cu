#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "../include/commit_gpu.cuh"
#include "../include/merkle.cuh"
#include "../include/hash.cuh"
#include "../include/field.cuh"
#define FIELD_WORDS 5
#define HASH_WORDS 4


__device__ void print_field_kernel(const char *label, const uint64_t *field, int field_words) {
    printf("%s: ", label);
    for (int i = 0; i < field_words; i++) {
        printf("%016lx ", field[i]);
    }
    printf("\n");
}

__host__ void print_field_host(const char *label, const uint64_t *field, int field_words) {
    printf("%s: ", label);
    for (int i = 0; i < field_words; i++) {
        printf("%016lx ", field[i]);
    }
    printf("\n");
}

// Kernel for computing layers with codewords up to layer 11
__global__ void commit_kernel(
    uint64_t *device_codeword, uint64_t *device_codeword_nxt, 
    uint64_t *device_alpha, uint64_t *device_offset,
    uint64_t *device_denominator_inv, uint64_t *device_eval_basis,
    uint64_t *device_temp1, uint64_t *device_temp2,
    uint64_t *device_temp3, uint64_t *device_temp4,
    uint64_t *device_temp5, uint64_t *device_alpha_offset, uint64_t *device_layer_hashes,
    uint64_t *device_concatenated_tree, int N, int basis_len, int layer
) {
    int I = blockIdx.x * blockDim.x + threadIdx.x;
    if (I >= N / 2) return;

    int idx1 = 2 * I * FIELD_WORDS;
    int idx2 = (2 * I + 1) * FIELD_WORDS;
    int idx3 = I * FIELD_WORDS;  // Used for codeword computation
    int tree_idx = I * (1 + FIELD_WORDS + HASH_WORDS);  // For [index, hash, element] in concatenated tree

    field_sub(&device_temp1[idx3], &device_codeword[idx1], &device_codeword[idx2], FIELD_WORDS);
    i_th_ele_in_span(&device_alpha_offset[idx3], device_eval_basis, basis_len, 2 * I);
    field_addEqual(&device_alpha_offset[idx3], device_offset, FIELD_WORDS);

    field_sub(&device_temp2[idx3], device_alpha, &device_alpha_offset[idx3], FIELD_WORDS);
    field_mul(&device_temp3[idx3], &device_temp2[idx3], device_denominator_inv, FIELD_WORDS);
    field_mul(&device_temp4[idx3], &device_temp3[idx3], &device_temp1[idx3], FIELD_WORDS);
    field_add(&device_temp5[idx3], &device_temp4[idx3], &device_codeword[idx1], FIELD_WORDS);

    memcpy(&device_codeword_nxt[idx3], &device_temp5[idx3], FIELD_WORDS * sizeof(uint64_t));

    uint64_t combined[FIELD_WORDS + HASH_WORDS];
    uint8_t digest[HASH_SIZE];

    if (layer == 0) {
        // Layer 0: initial hashes of codewords
        memcpy(combined, &device_codeword[idx1], FIELD_WORDS * sizeof(uint64_t));
        memcpy(combined + FIELD_WORDS, &device_codeword[idx2], FIELD_WORDS * sizeof(uint64_t));
        SHA3(digest, (uint8_t *)combined, 2 * FIELD_WORDS * sizeof(uint64_t), 256);
        memcpy(&device_concatenated_tree[tree_idx + 1], digest, HASH_WORDS * sizeof(uint64_t));
    } else if (layer > 0 && layer <= 11) {
        // Intermediate layers: combine codewords and previous layer hashes
        memcpy(combined, &device_codeword[idx1], FIELD_WORDS * sizeof(uint64_t));
        memcpy(combined + FIELD_WORDS, &device_concatenated_tree[tree_idx - HASH_WORDS], HASH_WORDS * sizeof(uint64_t));
        SHA3(digest, (uint8_t *)combined, (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t), 256);
        memcpy(&device_concatenated_tree[tree_idx], digest, HASH_WORDS * sizeof(uint64_t));
        memcpy(&device_concatenated_tree[tree_idx + HASH_WORDS], &device_codeword[idx3], FIELD_WORDS * sizeof(uint64_t));
    }
    __syncthreads();
}

// Kernel for computing remaining layers with only hashes until the Merkle root
__global__ void merkle_kernel(
    uint64_t *device_layer_hashes, 
    uint64_t *device_merkle_root, 
    uint64_t *device_concatenated_tree, 
    int N
) {
    int I = blockIdx.x * blockDim.x + threadIdx.x;

    while (N > 1) {
        if (I < N / 2) {
            uint64_t combined[2 * HASH_WORDS];
            uint8_t digest[HASH_SIZE];

            int idx1 = (2 * I) * HASH_WORDS;
            int idx2 = (2 * I + 1) * HASH_WORDS;
            int idx3 = I * HASH_WORDS;

            memcpy(combined, &device_layer_hashes[idx1], HASH_WORDS * sizeof(uint64_t));
            memcpy(combined + HASH_WORDS, &device_layer_hashes[idx2], HASH_WORDS * sizeof(uint64_t));
            SHA3(digest, (uint8_t *)combined, 2 * HASH_WORDS * sizeof(uint64_t), 256);
            memcpy(&device_layer_hashes[idx3], digest, HASH_WORDS * sizeof(uint64_t));

            int concat_index = I * (1 + HASH_WORDS);
            device_concatenated_tree[concat_index] = idx3;  // Storing index
            memcpy(&device_concatenated_tree[concat_index + 1], digest, HASH_WORDS * sizeof(uint64_t));
        }

        __syncthreads();
        N /= 2;
    }

    if (I == 0) {
        memcpy(device_merkle_root, device_layer_hashes, HASH_SIZE);
        printf("Final Merkle root: ");
        for (int j = 0; j < HASH_WORDS; ++j) {
            printf("%016lx ", device_merkle_root[j]);
        }
        printf("\n");
    }
}


void commit_launch(
    uint64_t **codeword, uint64_t **codeword_nxt, 
    uint64_t *alpha, uint64_t *offset, 
    uint64_t denominator_inv, uint64_t *eval_basis, 
    int N, uint64_t *root, merkleTree *tree 
) {
    printf("Starting commit_launch\n");
    printf("N = %d, FIELD_WORDS = %d\n", N, FIELD_WORDS);
    int basis_len = (int)log2(N);
    printf("basis len: %d\n", basis_len);

    int layer = 17 - (int)log2(N);
//     if(N == 131072) {
//     for (int i = 0; i < N; ++i) {
//         tree[i].index = i + 1;
//         memcpy(tree[i].element, codeword[i], FIELD_WORDS * sizeof(uint64_t));
//         memset(tree[i].hash, 0, HASH_WORDS * sizeof(uint64_t));
//         tree[i].tree_layer = 0;
//     }
// }

    // if (N == 131072) {
    //     initialize_file("temp1.txt");
    //     initialize_file("temp2.txt");
    //     initialize_file("temp3.txt");
    //     initialize_file("temp4.txt");
    //     initialize_file("temp5.txt");
    //     initialize_file("alpha_offset.txt");
    // }

    uint64_t *flattened_codeword = (uint64_t *)malloc(N * FIELD_WORDS * sizeof(uint64_t));
    if (flattened_codeword == NULL) {
        fprintf(stderr, "Error: malloc failed for flattened_codeword\n");
        return;
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < FIELD_WORDS; ++j) {
            int index = i * FIELD_WORDS + j;
            flattened_codeword[index] = codeword[i][j];
        }
    }

    printf("First few flattened codeword values:\n");
    for (int i = 0; i < 10; ++i) {
        printf("%016lx ", flattened_codeword[i]);
    }
    printf("\n");

    uint64_t *flattened_codeword_nxt = (uint64_t *)malloc((N / 2) * FIELD_WORDS * sizeof(uint64_t));
    if (flattened_codeword_nxt == NULL) {
        fprintf(stderr, "Error: malloc failed for flattened_codeword_nxt\n");
        free(flattened_codeword);
        return;
    }

    for (int i = 0; i < N / 2; ++i) {
        codeword_nxt[i] = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
        if (codeword_nxt[i] == NULL) {
            fprintf(stderr, "Error: malloc failed for codeword_nxt[%d]\n", i);
            for (int j = 0; j < i; ++j) {
                free(codeword_nxt[j]);
            }
            free(flattened_codeword);
            free(flattened_codeword_nxt);
            return;
        }
    }

    uint64_t flattened_eval_basis[basis_len];
    for (int i = 0; i < basis_len; ++i) {
        flattened_eval_basis[i] = eval_basis[i];
    }

    int field_size = N * FIELD_WORDS * sizeof(uint64_t);

    uint64_t *device_codeword, *device_codeword_nxt, *device_alpha, *device_offset;
    uint64_t *device_denominator_inv, *device_eval_basis;
    uint64_t *device_temp1, *device_temp2, *device_temp3, *device_temp4, *device_temp5, *device_alpha_offset;
    uint64_t *device_layer_hashes, *device_merkle_root, *device_concatenated_tree;
    uint64_t *flattened_alpha_offset = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp1 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp2 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp3 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp4 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp5 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));

    flattened_temp1[N/2 * FIELD_WORDS] = {0}, flattened_temp2[N/2 * FIELD_WORDS] = {0}, flattened_temp3[N/2 * FIELD_WORDS] = {0}, flattened_temp4[N/2 * FIELD_WORDS] = {0}, flattened_temp5[N/2 * FIELD_WORDS] = {0};
    flattened_alpha_offset[N/2 * FIELD_WORDS] = {0};
    cudaMalloc((void**)&device_codeword, field_size);
    cudaMalloc((void**)&device_codeword_nxt, field_size);
    cudaMalloc((void**)&device_alpha, FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void**)&device_offset, FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void**)&device_denominator_inv, sizeof(uint64_t));
    cudaMalloc((void**)&device_eval_basis, basis_len * sizeof(uint64_t));
    cudaMalloc((void **)&device_temp1, N/2 * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void **)&device_temp2, N/2 * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void **)&device_temp3, N/2 * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void **)&device_temp4, N/2 * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void **)&device_temp5, N/2 * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void **)&device_alpha_offset, N/2 * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void**)&device_layer_hashes, N * HASH_WORDS * sizeof(uint64_t));
    cudaMalloc((void**)&device_merkle_root, HASH_SIZE);
    cudaMalloc((void**)&device_concatenated_tree, (N / 2) * (1 + FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));

    cudaMemcpy(device_codeword, flattened_codeword, field_size, cudaMemcpyHostToDevice);
    cudaMemcpy(device_alpha, alpha, FIELD_WORDS * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(device_offset, offset, sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(device_denominator_inv, &denominator_inv, sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(device_eval_basis, flattened_eval_basis, basis_len * sizeof(uint64_t), cudaMemcpyHostToDevice);

    int threads_per_block = 1;
    int num_blocks = (N / 2 + threads_per_block - 1) / threads_per_block;
    commit_kernel<<<num_blocks * 2, threads_per_block>>>(
        device_codeword, device_codeword_nxt, device_alpha, device_offset,
        device_denominator_inv, device_eval_basis, device_temp1, device_temp2,
        device_temp3, device_temp4, device_temp5, device_alpha_offset, device_layer_hashes, device_concatenated_tree, N, basis_len, layer
    );
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Kernel launch failed: %s\n", cudaGetErrorString(err));
        free(flattened_codeword);
        free(flattened_codeword_nxt);
        free(flattened_eval_basis);
        return;
    }
    // Retrieve data from device_concatenated_tree and populate the `tree` structure
    uint64_t *flattened_tree = (uint64_t *)malloc((N / 2) * (1 + FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
    cudaMemcpy(flattened_tree, device_concatenated_tree, (N / 2) * (1 + FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_codeword_nxt, device_codeword_nxt, (N / 2) * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp1, device_temp1, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp2, device_temp2, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp3, device_temp3, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp4, device_temp4, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp5, device_temp5, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_alpha_offset, device_alpha_offset, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);

    // Populate the `tree` structure with unflattened data
    for (int i = 0; i < N / 2; ++i) {
        int offset = i * (1 + FIELD_WORDS + HASH_WORDS);
        tree[N / 2 + i].index = flattened_tree[offset];  // First element is the index
        memcpy(tree[N / 2 + i].element, &flattened_tree[offset + 1], FIELD_WORDS * sizeof(uint64_t));  // Next FIELD_WORDS elements are the element values
        memcpy(tree[N / 2 + i].hash, &flattened_tree[offset + 1 + FIELD_WORDS], HASH_WORDS * sizeof(uint64_t));  // Next HASH_WORDS elements are the hash
        tree[N / 2 + i].tree_layer = layer + 1;
    }

    // write_to_file("temp1.txt", flattened_temp1, FIELD_WORDS, N/2);
    // write_to_file("temp2.txt", flattened_temp2, FIELD_WORDS, N/2);
    // write_to_file("temp3.txt", flattened_temp3, FIELD_WORDS, N/2);
    // write_to_file("temp4.txt", flattened_temp4, FIELD_WORDS, N/2);
    // write_to_file("temp5.txt", flattened_temp5, FIELD_WORDS, N/2);
    // write_to_file("alpha_offset.txt", flattened_alpha_offset, FIELD_WORDS, N/2);

    if (N == 32) {
        while (N > 1) {
            printf("reaching here tanjan\n");
            int tpb = 32; 
            int nb = 1;
            merkle_kernel<<<nb, tpb>>>(device_layer_hashes, device_merkle_root, device_concatenated_tree, N);
            N = N / 2;
            cudaDeviceSynchronize();
            cudaMemcpy(root, device_merkle_root, HASH_SIZE, cudaMemcpyDeviceToHost);
            for (int i = 0; i < N / 2; ++i) {
                int offset = i * (1 + FIELD_WORDS + HASH_WORDS);
                tree[N / 2 + i].index = flattened_tree[offset];
                memcpy(tree[N / 2 + i].element, &flattened_tree[offset + 1], FIELD_WORDS * sizeof(uint64_t));
                memcpy(tree[N / 2 + i].hash, &flattened_tree[offset + 1 + FIELD_WORDS], HASH_WORDS * sizeof(uint64_t));
                tree[N / 2 + i].tree_layer = layer + 1;
            }
        }
    }

    for (int i = 0; i < N / 2; ++i) {
        for (int j = 0; j < FIELD_WORDS; ++j) {
            int index = i * FIELD_WORDS + j;
            if (index >= (N / 2) * FIELD_WORDS) {
                printf("Out-of-bounds access at index: %d\n", index);
                break;
            }
            codeword_nxt[i][j] = flattened_codeword_nxt[index];
        }
    }

    printf("First few codeword_nxt values:\n");
    for (int i = 0; i < 10; i++) {
        printf("%016lx ", flattened_codeword_nxt[i]);
    }
    printf("\n");

    cudaFree(device_codeword);
    cudaFree(device_codeword_nxt);
    cudaFree(device_alpha);
    cudaFree(device_offset);
    cudaFree(device_denominator_inv);
    cudaFree(device_eval_basis);
    cudaFree(device_merkle_root);
    cudaFree(device_layer_hashes);
    cudaFree(device_concatenated_tree);
    free(flattened_codeword);
    free(flattened_codeword_nxt);
    free(flattened_tree);

    printf("Memory freed and commit_launch completed.\n");
}