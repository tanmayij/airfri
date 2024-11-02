#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "../include/commit-launch-merkle.cuh"
#include "../include/hash-host.cuh"
#include "../include/field.cuh"
#include "../include/merkle.cuh"

__host__ __device__ void i_th_ele_in_span_commit(uint64_t *result, uint64_t *basis, int len_basis, int i)
{
    *result = 0;  
    for (int bit = 0; bit < len_basis; ++bit)
    {
        if ((i >> bit) & 1)
        {
            *result ^= basis[bit];  
        }
    }
}

const size_t field_words = 4;

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

__host__ void initialize_file(const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return;
    }
    fclose(file); 
}

__host__ void write_to_file(const char *filename, const uint64_t *data, int field_words, int total_indices) {
    FILE *file = fopen(filename, "a");
    if (file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return;
    }

    for (int index = 0; index < total_indices; index++) {
        fprintf(file, "Index %d: ", index);
        for (int i = 0; i < field_words; i++) {
            fprintf(file, "%llx ", (unsigned long long)data[index + i]);
        }
        fprintf(file, "\n");
    }
    fflush(file);  
    fclose(file);
}

__global__ void commit_kernel(
    uint64_t *device_codeword, uint64_t *device_codeword_nxt,
    uint64_t *device_alpha, uint64_t *device_offset,
    uint64_t *device_denominator_inv, uint64_t *device_eval_basis,
    uint64_t *device_temp1, uint64_t *device_temp2,
    uint64_t *device_temp3, uint64_t *device_temp4,
    uint64_t *device_temp5, uint64_t *device_alpha_offset, 
    uint64_t *device_layer_hashes, uint64_t *device_concatenated_tree, int N, int basis_len
) {
    size_t I = blockIdx.x * blockDim.x + threadIdx.x;

    if (I >= N / 2) return;

    int idx1 = 2 * I * FIELD_WORDS;
    int idx2 = (2 * I + 1) * FIELD_WORDS;
    int idx3 = I * FIELD_WORDS;

    field_sub(&device_temp1[idx3], &device_codeword[idx1], &device_codeword[idx2], field_words);
    i_th_ele_in_span_commit(&device_alpha_offset[idx3], device_eval_basis, basis_len, 2 * I);
    field_addEqual(&device_alpha_offset[idx3], device_offset, field_words); 

    field_sub(&device_temp2[idx3], device_alpha, &device_alpha_offset[idx3], field_words);
    field_mul(&device_temp3[idx3], &device_temp2[idx3], device_denominator_inv, field_words);
    field_mul(&device_temp4[idx3], &device_temp3[idx3], &device_temp1[idx3], field_words);
    field_add(&device_temp5[idx3], &device_temp4[idx3], &device_codeword[idx1], field_words);

    memcpy(&device_codeword_nxt[idx3], &device_temp5[idx3], field_words * sizeof(uint64_t));

    uint64_t combined[2 * HASH_WORDS];
    uint8_t digest[HASH_SIZE];
    if (I < N / 2 && N == 8) { //hardcoded for now
        memcpy(combined, &device_codeword[idx1], FIELD_WORDS * sizeof(uint64_t));  
        memcpy(combined + FIELD_WORDS, &device_codeword[idx2], FIELD_WORDS * sizeof(uint64_t));  
        SHA3(digest, (uint8_t *)combined, (FIELD_WORDS + FIELD_WORDS) * sizeof(uint64_t), 256);
        memcpy((uint8_t *)&device_layer_hashes[I * HASH_WORDS], digest, 32);
        memcpy(&device_concatenated_tree[idx3], combined, CONCAT_WORDS * sizeof(uint64_t));
    }
    __syncthreads();

    for (int layer = 1; layer < 3; layer++) {
        if (I < (N / (2 * (1 << layer)))) {
            int codeword_idx = (1 << layer) * idx3;
            memcpy(combined, &device_layer_hashes[(I * 2) * HASH_WORDS], HASH_WORDS * sizeof(uint64_t));
            memcpy(combined + FIELD_WORDS, &device_codeword[codeword_idx], FIELD_WORDS * sizeof(uint64_t));

            SHA3(digest, (uint8_t *)combined,  FIELD_WORDS+HASH_WORDS * sizeof(uint64_t), 256);
            memcpy((uint8_t *)&device_layer_hashes[I * HASH_WORDS], digest, 32);
            memcpy(&device_concatenated_tree[I * CONCAT_WORDS], combined, CONCAT_WORDS * sizeof(uint64_t));
        }
        __syncthreads();
    }
}

__global__ void merkle_kernel(
    uint64_t *device_layer_hashes, 
    uint64_t *device_merkle_root, 
    uint64_t *device_concatenated_tree, 
    int N
) {
    int I = blockIdx.x * blockDim.x + threadIdx.x;
    int layer = 0;

    // Traverse the tree up to the root, storing concatenated hashes in `device_concatenated_tree`
    while (N > 1) {
        if (I < N / 2) {
            uint64_t combined[2 * HASH_WORDS];
            uint8_t digest[HASH_SIZE];

            // Combine two hashes from the current layer
            memcpy(combined, &device_layer_hashes[(2 * I) * HASH_WORDS], HASH_WORDS * sizeof(uint64_t));
            memcpy(combined + HASH_WORDS, &device_layer_hashes[(2 * I + 1) * HASH_WORDS], HASH_WORDS * sizeof(uint64_t));

            // Hash the combined result to obtain the parent node
            SHA3(digest, (uint8_t *)combined, 2 * HASH_WORDS * sizeof(uint64_t), 256);

            // Store the digest in the parent position in the next layer
            memcpy((uint8_t *)&device_layer_hashes[I * HASH_WORDS], digest, 32);

            // Store the combined hash values in `device_concatenated_tree` for this layer
            int concat_index = (layer * N / 2 + I) * CONCAT_WORDS;
            memcpy(&device_concatenated_tree[concat_index], combined, 2 * HASH_WORDS * sizeof(uint64_t));
        }

        // Synchronize threads before moving to the next layer and halve the number of nodes
        __syncthreads();
        N /= 2;
        layer++;
    }

    // Assign the root hash to the output only once, when N == 1
    if (I == 0) {
        memcpy(device_merkle_root, device_layer_hashes, HASH_SIZE);
        printf("Final Merkle root: ");
        for (int j = 0; j < HASH_WORDS; ++j) {
            printf("%016lx ", device_merkle_root[j]);
        }
        printf("\n");
    }
}
// __global__ void merkle_kernel(uint64_t *device_layer_hashes, uint64_t *device_merkle_root, int N) {
//     int I = blockIdx.x * blockDim.x + threadIdx.x;

//     while (N > 1) {
//         if (I < N / 2) {
//             uint64_t combined[2 * HASH_WORDS];
//             uint8_t digest[HASH_SIZE];

//             memcpy(combined, &device_layer_hashes[(2 * I) * HASH_WORDS], HASH_WORDS * sizeof(uint64_t));  
//             memcpy(combined + HASH_WORDS, &device_layer_hashes[(2 * I + 1) * HASH_WORDS], HASH_WORDS * sizeof(uint64_t));  

//             SHA3(digest, (uint8_t *)combined, 2 * HASH_WORDS * sizeof(uint64_t), 256);
//             memcpy((uint8_t *)&device_layer_hashes[I * HASH_WORDS], digest, 32);
//         }

//         __syncthreads();
//         N /= 2;
//     }

//     if (I == 0) {
//         memcpy(device_merkle_root, device_layer_hashes, HASH_SIZE);  
//         printf("Final Merkle root: ");
//         for (int j = 0; j < HASH_WORDS; ++j) {
//             printf("%016lx ", device_merkle_root[j]);
//         }
//         printf("\n");
//     }
// }

void commit_launch(
    uint64_t **codeword, uint64_t **codeword_nxt, 
    uint64_t *alpha, uint64_t *offset, 
    uint64_t denominator_inv, uint64_t *eval_basis, 
    int N, uint64_t *root, uint64_t *tree
) {
    printf("Starting commit_launch for toy example\n");
    int basis_len = (int)log2(N);

    uint64_t *flattened_codeword = (uint64_t *)malloc(N * FIELD_WORDS * sizeof(uint64_t));
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < FIELD_WORDS; ++j) {
            flattened_codeword[i * FIELD_WORDS + j] = codeword[i][j];
        }
    }

    uint64_t *device_codeword, *device_codeword_nxt, *device_alpha, *device_offset;
    uint64_t *device_denominator_inv, *device_eval_basis;
    uint64_t *device_temp1, *device_temp2, *device_temp3, *device_temp4, *device_temp5, *device_alpha_offset;
    uint64_t *device_layer_hashes, *device_merkle_root, *device_concatenated_tree;

    cudaMalloc((void**)&device_codeword, N * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void**)&device_codeword_nxt, (N / 2) * FIELD_WORDS * sizeof(uint64_t));
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
    cudaMalloc((void**)&device_concatenated_tree, (2 * N - 1) * CONCAT_WORDS * sizeof(uint64_t));

    cudaMemcpy(device_codeword, flattened_codeword, N * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(device_alpha, alpha, FIELD_WORDS * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(device_offset, offset, FIELD_WORDS * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(device_denominator_inv, &denominator_inv, sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(device_eval_basis, eval_basis, basis_len * sizeof(uint64_t), cudaMemcpyHostToDevice);

    int threads_per_block = 1;
    int num_blocks = (N / 2 + threads_per_block - 1) / threads_per_block;
    commit_kernel<<<num_blocks, threads_per_block>>>(
        device_codeword, device_codeword_nxt, device_alpha, device_offset,
        device_denominator_inv, device_eval_basis, device_temp1, device_temp2,
        device_temp3, device_temp4, device_temp5, device_alpha_offset, device_layer_hashes, device_concatenated_tree, N, basis_len
    );
    cudaDeviceSynchronize();

    while (N > 1) {
        int tpb = 8; 
        int nb = 1;
        merkle_kernel<<<nb, tpb>>>(device_layer_hashes, device_merkle_root, device_concatenated_tree, N);
        N = N / 2;
        cudaDeviceSynchronize();
        cudaMemcpy(root, device_merkle_root, HASH_SIZE, cudaMemcpyDeviceToHost);
        cudaMemcpy(tree, device_concatenated_tree, 4/*tree_size*/ * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    }

    free(flattened_codeword);
    printf("Commit launch and Merkle root computation complete.\n");
}