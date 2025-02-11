
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "../include/hash.cuh"
#include "../include/hash-host.cuh"
#include "../include/commit-launch-merkle.cuh"
#include "../include/field.cuh"
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
    uint64_t *device_layer_hashes, uint64_t *device_tree_layer, uint64_t *device_tree_layer_nxt, int N, int basis_len
) {
    size_t I = blockIdx.x * blockDim.x + threadIdx.x;

    if (I >= N / 2) return;

    int idx1 = 2 * I * FIELD_WORDS;
    int idx2 = (2 * I + 1) * FIELD_WORDS;
    int idx3 = I * FIELD_WORDS;

    field_sub(&device_temp1[idx3], &device_codeword[idx1], &device_codeword[idx2], field_words);
    i_th_ele_in_span(&device_alpha_offset[idx3], device_eval_basis, basis_len, 2 * I);
    field_addEqual(&device_alpha_offset[idx3], device_offset, field_words); 

    field_sub(&device_temp2[idx3], device_alpha, &device_alpha_offset[idx3], field_words);
    field_mul(&device_temp3[idx3], &device_temp2[idx3], device_denominator_inv, field_words);
    field_mul(&device_temp4[idx3], &device_temp3[idx3], &device_temp1[idx3], field_words);
    field_add(&device_temp5[idx3], &device_temp4[idx3], &device_codeword[idx1], field_words);

    memcpy(&device_codeword_nxt[idx3], &device_temp5[idx3], field_words * sizeof(uint64_t));
}

__global__ void compute_tree_layers(uint64_t *device_codeword_nxt, uint64_t *device_layer_hashes, uint64_t *device_tree_layer,
uint64_t *device_tree_layer_nxt, uint64_t *device_combined_sibling_codewords, uint64_t *device_concat_codeword_to_hash, uint64_t *device_digest, int N)
{   
    size_t I = blockIdx.x * blockDim.x + threadIdx.x;

    if (I >= N / 2) return;
    if (I < N / 2 && N == 131072) {
        int idx1 = 2 * I * FIELD_WORDS;
        int idx2 = (2 * I + 1) * FIELD_WORDS;
        int idx3 = I * FIELD_WORDS;//for codeword_nxt element
        int idx5 = I * HASH_WORDS; //for hash
        int idx4 = I * (FIELD_WORDS + HASH_WORDS);
        //step 1
        memcpy(&device_combined_sibling_codewords[idx3], &device_tree_layer[idx1], FIELD_WORDS * sizeof(uint64_t));
        memcpy(&device_combined_sibling_codewords[idx3 + FIELD_WORDS], &device_tree_layer[idx2], FIELD_WORDS * sizeof(uint64_t));
        //step 2
        SHA3((uint8_t *)&device_digest[idx5], (uint8_t *)&device_combined_sibling_codewords[idx3], 2 * FIELD_WORDS * sizeof(uint64_t), 256);
        //step 3
        memcpy(&device_concat_codeword_to_hash[idx4], &device_codeword_nxt[idx3], FIELD_WORDS * sizeof(uint64_t));
        memcpy(&device_concat_codeword_to_hash[idx4 + FIELD_WORDS], &device_digest[idx5], HASH_WORDS * sizeof(uint64_t));
        //step 4: future-> remove this dependency
        memcpy(&device_tree_layer_nxt[idx4], &device_concat_codeword_to_hash[idx4], (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
    
    }

    if(I < N/2 && N < 131072 && N >= 64) {
        int idx1 = 2 * I * (FIELD_WORDS + HASH_WORDS);
        int idx2 = (2 * I + 1) * (FIELD_WORDS + HASH_WORDS);
        int idx3 = I * FIELD_WORDS;//for codeword_nxt element
        int idx5 = I * HASH_WORDS; //for hash
        int idx4 = I * (FIELD_WORDS + HASH_WORDS);
        int idx6 = I * 2 * (FIELD_WORDS + HASH_WORDS);
        //step 1
        memcpy(&device_combined_sibling_codewords[idx6], &device_tree_layer[idx1], (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
        memcpy(&device_combined_sibling_codewords[idx6 + (FIELD_WORDS + HASH_WORDS)], &device_tree_layer[idx2], (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
        //step 2
        SHA3((uint8_t *)&device_digest[idx5], (uint8_t *)&device_combined_sibling_codewords[idx6], 2 * (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t), 256);
        //step3
        memcpy(&device_concat_codeword_to_hash[idx4], &device_codeword_nxt[idx3], FIELD_WORDS * sizeof(uint64_t));
        memcpy(&device_concat_codeword_to_hash[idx4 + FIELD_WORDS], &device_digest[idx5], HASH_WORDS * sizeof(uint64_t));
        //step 4
        memcpy(&device_tree_layer_nxt[idx4], &device_concat_codeword_to_hash[idx4], (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));

    }
}
__global__ void merkle_kernel(
    uint64_t *device_layer_hashes, 
    uint64_t *device_merkle_root, 
    uint64_t *device_tree_layer,
    uint64_t *device_tree_layer_nxt,
    uint64_t *device_combined_sibling_codewords,
    uint64_t *device_digest,
    uint64_t *device_combined_sibling_hashes,
    int N
) {
    int I = blockIdx.x * blockDim.x + threadIdx.x;
    int idx1, idx2, idx3, idx4;

    if(I < N/2 && N == 32){ //only if N is 32, we use codewords
        idx1 = 2 * I * (FIELD_WORDS + HASH_WORDS);
        idx2 = (2 * I + 1) * (FIELD_WORDS + HASH_WORDS);
        idx3 = I * HASH_WORDS;
        idx4 = I * 2 * (FIELD_WORDS + HASH_WORDS);
        //step 1
        memcpy(&device_combined_sibling_codewords[idx4], &device_tree_layer[idx1], (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
        memcpy(&device_combined_sibling_codewords[idx4 + (FIELD_WORDS + HASH_WORDS)], &device_tree_layer[idx2], (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
        //step 2
        SHA3((uint8_t *)&device_digest[idx3],(uint8_t *)&device_combined_sibling_codewords[idx4], 2 * (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t), 256);
        //step 3 - save all hashes in tree_layer_nxt
        memcpy(&device_tree_layer_nxt[idx3], &device_digest[idx3], HASH_WORDS * sizeof(uint64_t) );
    }

    if(I < N/2 && N < 32) 
    {
        idx1 = 2 * I * (HASH_WORDS);
        idx2 = (2 * I + 1) * (HASH_WORDS);
        idx3 = I * HASH_WORDS;
        idx4 = I * 2 * HASH_WORDS;
        //step 1
        memcpy(&device_combined_sibling_hashes[idx4], &device_tree_layer[idx1], (HASH_WORDS) * sizeof(uint64_t));
        memcpy(&device_combined_sibling_hashes[idx4 + (HASH_WORDS)], &device_tree_layer[idx2], (HASH_WORDS) * sizeof(uint64_t));
        //step 2
        SHA3((uint8_t *)&device_digest[idx3], (uint8_t *)&device_combined_sibling_hashes[idx3], 2 * (HASH_WORDS) * sizeof(uint64_t), 256);
        //step 3 - save all hashes in tree_layer_nxt
        memcpy(&device_tree_layer_nxt[idx3], &device_digest[idx3], HASH_WORDS * sizeof(uint64_t) );

    }
    // if (I == 0 && N == 2) {
    //     // Combine the last two sibling hashes into device_combined_sibling_hashes
    //     int idx1 = 0 * HASH_WORDS;            // First sibling hash
    //     int idx2 = 1 * HASH_WORDS;            // Second sibling hash
    //     int idx3 = 0 * HASH_WORDS;            // Output digest index
    //     int idx4 = 0 * 2 * HASH_WORDS;        // Combined hash index
    
    //     // Copy the two sibling hashes into a combined buffer
    //     memcpy(&device_combined_sibling_hashes[idx4], &device_tree_layer[idx1], HASH_WORDS * sizeof(uint64_t));
    //     memcpy(&device_combined_sibling_hashes[idx4 + HASH_WORDS], &device_tree_layer[idx2], HASH_WORDS * sizeof(uint64_t));
    
    //     // Compute the SHA3 hash of the two combined sibling hashes
    //     SHA3((uint8_t *)&device_digest[idx3], (uint8_t *)&device_combined_sibling_hashes[idx4], 2 * HASH_WORDS * sizeof(uint64_t), 256);
    
    //     // Store the computed root in device_merkle_root
    //     memcpy(device_merkle_root, &device_digest[idx3], HASH_WORDS * sizeof(uint64_t));
    // }
    if (I == 0 && N == 2) {
        // Compute Merkle root within this kernel instead of launching another one
        int idx1 = 0 * HASH_WORDS; 
        int idx2 = 1 * HASH_WORDS;
        int idx3 = 0 * HASH_WORDS;
        int idx4 = 0 * 2 * HASH_WORDS;

        uint64_t combined_sibling_hashes[2 * HASH_WORDS];

        memcpy(&combined_sibling_hashes[0], &device_tree_layer[idx1], HASH_WORDS * sizeof(uint64_t));
        memcpy(&combined_sibling_hashes[HASH_WORDS], &device_tree_layer[idx2], HASH_WORDS * sizeof(uint64_t));

        SHA3((uint8_t *)&device_merkle_root[idx3], (uint8_t *)combined_sibling_hashes, 2 * HASH_WORDS * sizeof(uint64_t), 256);

        printf("Computed Merkle Root inside kernel: ");
        for (int j = 0; j < HASH_WORDS; j++) {
            printf("%016lx ", device_merkle_root[j]);
        }
        printf("\n");
    }
}

__global__ void compute_merkle_root_kernel(
    uint64_t *device_tree_layer,    // Input: layer with two sibling hashes
    uint64_t *device_merkle_root   // Output: the Merkle root
) {
    if (threadIdx.x == 0) {
        // Indices for the two sibling hashes
        int idx1 = 0 * HASH_WORDS;
        int idx2 = 1 * HASH_WORDS;

        // Combined buffer for the two sibling hashes
        uint64_t combined_sibling_hashes[2 * HASH_WORDS];

        // Copy the sibling hashes into the combined buffer
        memcpy(&combined_sibling_hashes[0], &device_tree_layer[idx1], HASH_WORDS * sizeof(uint64_t));
        memcpy(&combined_sibling_hashes[HASH_WORDS], &device_tree_layer[idx2], HASH_WORDS * sizeof(uint64_t));

        // Compute the Merkle root using SHA3
        SHA3(
            (uint8_t *)device_merkle_root,                  // Destination: Merkle root
            (uint8_t *)combined_sibling_hashes,            // Source: Combined sibling hashes
            2 * HASH_WORDS * sizeof(uint64_t),             // Input size: 2 sibling hashes
            256                                            // Output size: 256 bits
        );
    }
}

void commit_launch(
    uint64_t **codeword, uint64_t **codeword_nxt, 
    uint64_t *alpha, uint64_t *offset, 
    uint64_t denominator_inv, uint64_t *eval_basis, 
    int N, uint64_t *root, uint64_t **tree_layer, uint64_t **tree_layer_nxt, uint64_t ***tree
) {
    printf("Starting commit_launch\n");
    printf("N = %d, FIELD_WORDS = %d\n", N, FIELD_WORDS);
    int basis_len = (int)log2(N);
    printf("basis len: %d\n", basis_len);

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

    uint64_t *flattened_tree_layer;
    if(N == 131072) {
    flattened_tree_layer = (uint64_t *)malloc(N * FIELD_WORDS * sizeof(uint64_t));
    if (flattened_tree_layer == NULL) {
        fprintf(stderr, "Error: malloc failed for flattened_codeword\n");
        return;
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < FIELD_WORDS; ++j) {
            int index = i * FIELD_WORDS + j;
            flattened_tree_layer[index] = tree_layer[i][j];
        }
    }
    printf("First few flattened tree_layer values:\n");
    for (int i = 0; i < 10; ++i) {
        printf("%016lx ", flattened_tree_layer[i]);
    }
    printf("\n");

    } else {
    flattened_tree_layer = (uint64_t *)malloc(N * CONCAT_WORDS * sizeof(uint64_t));
    if (flattened_tree_layer == NULL) {
        fprintf(stderr, "Error: malloc failed for flattened_codeword\n");
        return;
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < CONCAT_WORDS; ++j) {
            int index = i * CONCAT_WORDS + j;
            flattened_tree_layer[index] = tree_layer[i][j];
        }
    }
    printf("First few flattened tree_layer values:\n");
    for (int i = 0; i < 10; ++i) {
        printf("%016lx ", flattened_tree_layer[i]);
    }
    printf("\n");
}


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

    uint64_t *flattened_tree_layer_nxt = (uint64_t *)malloc((N / 2) * CONCAT_WORDS * sizeof(uint64_t));
    if (flattened_tree_layer_nxt == NULL) {
        fprintf(stderr, "Error: malloc failed for flattened_tree_layer_nxt\n");
        free(flattened_codeword);
        return;
    }

    for (int i = 0; i < N / 2; ++i) {
        tree_layer_nxt[i] = (uint64_t *)malloc(CONCAT_WORDS * sizeof(uint64_t));
        if (tree_layer_nxt[i] == NULL) {
            fprintf(stderr, "Error: malloc failed for tree_layer_nxt[%d]\n", i);
            for (int j = 0; j < i; ++j) {
                free(tree_layer_nxt[j]);
            }
            free(flattened_tree_layer_nxt);
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
    uint64_t *device_layer_hashes, *device_merkle_root, *device_tree_layer, *device_tree_layer_nxt, *device_combined_sibling_codewords, *device_combined_sibling_hashes, *device_concat_codeword_to_hash, *device_digest;
    uint64_t *flattened_alpha_offset = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp1 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp2 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp3 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp4 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp5 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));

    flattened_temp1[N/2 * FIELD_WORDS] = {0}, flattened_temp2[N/2 * FIELD_WORDS] = {0}, flattened_temp3[N/2 * FIELD_WORDS] = {0}, flattened_temp4[N/2 * FIELD_WORDS] = {0}, flattened_temp5[N/2 * FIELD_WORDS] = {0};
    flattened_alpha_offset[N/2 * FIELD_WORDS] = {0};
    cudaMalloc((void**)&device_codeword, field_size);
    cudaMalloc((void**)&device_codeword_nxt, (N/2) * FIELD_WORDS * sizeof(uint64_t));
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
    if(N == 131072){
    cudaMalloc((void**)&device_tree_layer, field_size);
    cudaMalloc((void**)&device_combined_sibling_codewords, (N/2) * 2 * FIELD_WORDS * sizeof(uint64_t));
    } else 
    {
        cudaMalloc((void**)&device_tree_layer, N * CONCAT_WORDS * sizeof(uint64_t));
        cudaMalloc((void**)&device_combined_sibling_codewords, (N/2) * 2 * ( FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
    }
    cudaMalloc((void**)&device_tree_layer_nxt, (N/2) * CONCAT_WORDS * sizeof(uint64_t));
    cudaMalloc((void**)&device_merkle_root, HASH_WORDS * sizeof(uint64_t));
    cudaMalloc((void**)&device_concat_codeword_to_hash, (N/2) * (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
    cudaMalloc((void**)&device_combined_sibling_hashes, (N/2) * 2 * HASH_WORDS * sizeof(uint64_t));
    cudaMalloc((void**)&device_digest, (N/2) * HASH_WORDS * sizeof(uint64_t));

    cudaMemcpy(device_codeword, flattened_codeword, field_size, cudaMemcpyHostToDevice);
    cudaMemcpy(device_alpha, alpha, FIELD_WORDS * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(device_offset, offset, sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(device_denominator_inv, &denominator_inv, sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(device_eval_basis, flattened_eval_basis, basis_len * sizeof(uint64_t), cudaMemcpyHostToDevice);
    if(N == 131072) {
    cudaMemcpy(device_tree_layer, flattened_tree_layer,  N * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyHostToDevice);
    } else 
    {
        cudaMemcpy(device_tree_layer, flattened_tree_layer,  N * CONCAT_WORDS * sizeof(uint64_t), cudaMemcpyHostToDevice);
    }
    int threads_per_block = 1;
    int num_blocks = (N / 2 + threads_per_block - 1) / threads_per_block;
    commit_kernel<<<num_blocks * 2, threads_per_block>>>(
        device_codeword, device_codeword_nxt, device_alpha, device_offset,
        device_denominator_inv, device_eval_basis, device_temp1, device_temp2,
        device_temp3, device_temp4, device_temp5, device_alpha_offset, device_layer_hashes, device_tree_layer, device_tree_layer_nxt, N, basis_len
    );
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Kernel launch failed: %s\n", cudaGetErrorString(err));
        free(flattened_codeword);
        free(flattened_codeword_nxt);
        // free(flattened_tree_layer);
        //free(flattened_eval_basis);
        return;
    }
    
    compute_tree_layers<<<num_blocks * 2, threads_per_block>>> (
        device_codeword_nxt, device_layer_hashes, device_tree_layer, device_tree_layer_nxt, device_combined_sibling_codewords, device_concat_codeword_to_hash, device_digest, N
    );
    cudaDeviceSynchronize();

    err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Kernel launch failed: %s\n", cudaGetErrorString(err));
        free(flattened_codeword);
        free(flattened_codeword_nxt);
        // free(flattened_tree_layer);
        //free(flattened_eval_basis);
        return;
    }

    cudaMemcpy(flattened_codeword_nxt, device_codeword_nxt, (N / 2) * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp1, device_temp1, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp2, device_temp2, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp3, device_temp3, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp4, device_temp4, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp5, device_temp5, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_alpha_offset, device_alpha_offset, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_tree_layer_nxt, device_tree_layer_nxt, (N / 2) * CONCAT_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);

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

    for (int i = 0; i < N / 2; ++i) {
        for (int j = 0; j < CONCAT_WORDS; ++j) {
            int index = i * CONCAT_WORDS + j;
            if (index >= (N / 2) * CONCAT_WORDS) {
                printf("Out-of-bounds access at index: %d\n", index);
                break;
            }
            tree_layer_nxt[i][j] = flattened_tree_layer_nxt[index];
        }
    }

    printf("First few tree_layer_nxt values:\n");
    for (int i = 0; i < 10; i++) {
        printf("%016lx ", flattened_tree_layer_nxt[i]);
    }
    printf("\n");
    cudaFree(device_tree_layer);
    cudaFree(device_tree_layer_nxt);
    free(flattened_tree_layer_nxt);
    // write_to_file("temp1.txt", flattened_temp1, FIELD_WORDS, N/2);
    // write_to_file("temp2.txt", flattened_temp2, FIELD_WORDS, N/2);
    // write_to_file("temp3.txt", flattened_temp3, FIELD_WORDS, N/2);
    // write_to_file("temp4.txt", flattened_temp4, FIELD_WORDS, N/2);
    // write_to_file("temp5.txt", flattened_temp5, FIELD_WORDS, N/2);
    // write_to_file("alpha_offset.txt", flattened_alpha_offset, FIELD_WORDS, N/2);

    // if (N == 64) { 
    //     int tree_idx = 2 * (int)log2(N);
    //     tree[tree_idx] = (uint64_t **)malloc((N/2) * sizeof(uint64_t *));
    //     for(int i = 0; i < N/2; i ++){
    //         tree[tree_idx][i] = (uint64_t *)malloc(CONCAT_WORDS * sizeof(uint64_t));
    //     }
    //     tree_idx++; 
    //     int next_N = N / 2;
    //     for(int layer = tree_idx; layer <= 18; layer++){
    //         tree[layer] = (uint64_t **)malloc((next_N / 2) * sizeof(uint64_t *));
    //         for(int i = 0; i < next_N / 2; i++){
    //             tree[layer][i] = (uint64_t *)malloc(2 * HASH_WORDS * sizeof(uint64_t));
    //             if (!tree[layer][i]) {
    //                 printf("Memory allocation failed at tree[%d][%d]\n", layer, i);
    //                 exit(1); 
    //             }
    //         }
    //         next_N = next_N / 2;
    //         printf("Allocated memory for tree[%d] with %d elements\n", layer, (next_N / 2) * HASH_WORDS);
    //     }
    
    //     int index = 13;
    //     int threads_per_block = 32;
    //     int num_blocks = (N / 2 + threads_per_block - 1) / threads_per_block;
    
    //     // Copy the last flattened layer from device to host
    //     cudaMemcpy(flattened_tree_layer_nxt, device_tree_layer_nxt, (N / 2) * HASH_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    
    //     // Unflatten and assign to tree[12]
    //     for (int i = 0; i < N / 2; i++) {
    //         for (int j = 0; j < CONCAT_WORDS; j++) {
    //             tree[12][i][j] = flattened_tree_layer_nxt[i * CONCAT_WORDS + j];
    //         }
    //     }
    
    //     free(flattened_tree_layer_nxt);
    //     cudaFree(device_tree_layer_nxt);
    
    //     printf("Starting iterative Merkle kernel computation up the tree\n");
    //     N /= 2;
    
    //     // Allocate for iterative usage
    //     flattened_tree_layer_nxt = (uint64_t *)malloc((N / 2) * CONCAT_WORDS * sizeof(uint64_t));
    //     cudaMalloc((void **)&device_tree_layer_nxt, (N / 2) * CONCAT_WORDS * sizeof(uint64_t));
    
    //     while (N > 1) {
    //         int tpb = N / 2;
    //         int nb = (N + tpb - 1) / tpb;
    
    //         merkle_kernel<<<nb, tpb>>>(device_layer_hashes, device_merkle_root, device_tree_layer, device_tree_layer_nxt, N);
    //         cudaDeviceSynchronize();
    
    //         // Copy next layer from device to host
    //         cudaMemcpy(flattened_tree_layer_nxt, device_tree_layer_nxt, (N / 2) * HASH_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    
    //         for (int i = 0; i < N / 2; i++) {
    //             for (int j = 0; j < HASH_WORDS; j++) {
    //                 tree[index][i][j] = flattened_tree_layer_nxt[i * HASH_WORDS + j];
    //             }
    //         }
    
    //         if (N == 2) {
    //             cudaMemcpy(root, device_merkle_root, HASH_SIZE, cudaMemcpyDeviceToHost);
    //         }
    
    //         N /= 2;
    //         index++;
    //     }
    
    //     free(flattened_tree_layer_nxt);
    //     cudaFree(device_tree_layer_nxt);
    // }
    //something's up here
    if (N == 64) { 
        int tree_idx = 12;  //start with layer 12 for N == 64
        int next_N = N / 2; //initialize to 32 for the next layer size
        
        //allocate flattened memory for device -> host transfer
        uint64_t *flattened_tree_layer_nxt = (uint64_t *)malloc((N / 2) * CONCAT_WORDS * sizeof(uint64_t));
    
        //step 1: unflatten tree_layer_nxt computed by commit_kernel and assign it to tree[12]
        cudaMalloc((void **)&device_tree_layer_nxt, (N / 2) * CONCAT_WORDS * sizeof(uint64_t));
        cudaMemcpy(flattened_tree_layer_nxt, device_tree_layer_nxt, (N / 2) * CONCAT_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    
        tree[tree_idx] = (uint64_t **)malloc((next_N) * sizeof(uint64_t *));
        for (int i = 0; i < next_N; i++) {
            tree[tree_idx][i] = (uint64_t *)malloc(CONCAT_WORDS * sizeof(uint64_t));
            for (int j = 0; j < CONCAT_WORDS; j++) {
                tree[tree_idx][i][j] = flattened_tree_layer_nxt[i * CONCAT_WORDS + j];
            }
        }
        tree_idx++;  // Move to the next tree layer index

        // Step 2: Assign tree_layer_nxt to tree_layer for the upcoming Merkle kernel computation
        cudaMemcpy(device_tree_layer, device_tree_layer_nxt, (next_N / 2) * CONCAT_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToDevice);
        // N = N/2;//32
        // //step3: loop over the remaining layers and only compute hashes.
        // while(N>=2){ //reduce until N becomes 2
        //     int tpb = N;
        //     if(tpb == 0) {tpb = 1;}
        //     int nb = (N + (tpb - 1)) / tpb;

        //     //launch merkle_kernel which will prallely compute the merkle tree layers for each layer. start when layer has 32 elements
        //     merkle_kernel<<<nb, tpb>>>(device_layer_hashes, device_merkle_root, device_tree_layer, device_tree_layer_nxt, device_combined_sibling_codewords, device_digest, device_combined_sibling_hashes, N);
        //     cudaDeviceSynchronize();

        //     //copy the next layer in the tree from device to host
        //     cudaMalloc((void **)&device_tree_layer_nxt, (N / 2) * HASH_WORDS * sizeof(uint64_t)); //next layer in devicewould have N/2 elements
        //     cudaMemcpy(flattened_tree_layer_nxt, device_tree_layer_nxt, (N/2) * HASH_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToDevice); //copy to device_layer_nxt
        //     printf("First few tree_layer_nxt values in the looooop:\n");
        //     for (int i = 0; i < HASH_WORDS; i++) {
        //         printf("%016lx ", flattened_tree_layer_nxt[i]);
        //     }
        //     printf("\n");
        //     //unflatten and store tree layer values in 'tree'
        //     tree[tree_idx] = (uint64_t **)malloc((next_N / 2) * sizeof(uint64_t *));
        //     for (int i = 0; i < next_N / 2; i++) {
        //         tree[tree_idx][i] = (uint64_t *)malloc(HASH_WORDS * sizeof(uint64_t));
        //         for (int j = 0; j < HASH_WORDS; j++) {
        //             tree[tree_idx][i][j] = flattened_tree_layer_nxt[i * HASH_WORDS + j];
        //         }
        //     }
        //     //update tree_idx? 
        // }
    
        //step 3: Loop over remaining layers, updating tree[layer] with each iteration
        while (next_N >= 2) {
            int tpb = min(32, next_N / 2);
            if(tpb == 0) {tpb = 1;}
            int nb = (next_N + tpb - 1) / tpb;
    
            // merkle_kernel for each layer (computes next layer hashes)
            merkle_kernel<<<nb, tpb>>>(device_layer_hashes, device_merkle_root, device_tree_layer, device_tree_layer_nxt, device_combined_sibling_codewords, device_digest, device_combined_sibling_hashes, next_N);
            cudaDeviceSynchronize();
            
            cudaMalloc((void **)&device_tree_layer_nxt, (next_N / 2) * HASH_WORDS * sizeof(uint64_t));
            cudaMemcpy(flattened_tree_layer_nxt, device_tree_layer_nxt, (next_N / 2) * HASH_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
            printf("First few tree_layer_nxt values in the looooop:\n");
            for (int i = 0; i < HASH_WORDS + 4; i++) {
                printf("%016lx ", flattened_tree_layer_nxt[i]);
            }
            printf("\n");
            //unflatten and store in tree[tree_idx]
            tree[tree_idx] = (uint64_t **)malloc((next_N / 2) * sizeof(uint64_t *));
            for (int i = 0; i < next_N / 2; i++) {
                tree[tree_idx][i] = (uint64_t *)malloc(HASH_WORDS * sizeof(uint64_t));
                for (int j = 0; j < HASH_WORDS; j++) {
                    tree[tree_idx][i][j] = flattened_tree_layer_nxt[i * HASH_WORDS + j];
                }
            }
            printf("Populated tree[%d] with %d elements\n", tree_idx, next_N / 2);
            tree_idx++;
            
            //step 4: Update device_tree_layer with the contents of device_tree_layer_nxt
            cudaMemcpy(device_tree_layer, device_tree_layer_nxt, (next_N / 2) * HASH_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToDevice);

            next_N /= 2;
        // Handle the final layer (tree[16]) to compute the Merkle root
        if (next_N == 2) {
            // Copy Merkle root directly from device after merkle_kernel execution
            cudaMemcpy(root, device_merkle_root, HASH_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
        
            // Print the computed Merkle root
            printf("Computed Merkle Root: ");
            for (int i = 0; i < HASH_WORDS; i++) {
                printf("%016lx ", root[i]);
            }
            printf("\n");
        }
        // if (next_N == 2) {
        //     // Allocate memory for the Merkle root on the device
        //     uint64_t *device_merkle_root;
        //     cudaMalloc((void **)&device_merkle_root, HASH_WORDS * sizeof(uint64_t));

        //     // Launch kernel to compute the Merkle root
        //     compute_merkle_root_kernel<<<1, 1>>>(
        //         device_tree_layer,     // Input: Final layer with two sibling hashes
        //         device_merkle_root     // Output: Merkle root
        //     );
        //     cudaDeviceSynchronize();

        //     // Copy the Merkle root from device to host
        //     cudaMemcpy(root, device_merkle_root, HASH_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);

        //     // Free device memory for the root
        //     cudaFree(device_merkle_root);

        //     // Print the computed Merkle root
        //     printf("Computed Merkle Root: ");
        //     for (int i = 0; i < HASH_WORDS; i++) {
        //         printf("%016lx ", root[i]);
        //     }
        //     printf("\n");
        // }
        } //close while loop
        free(flattened_tree_layer_nxt);
    }

    
    
    cudaFree(device_codeword_nxt);
    cudaFree(device_alpha);
    cudaFree(device_offset);
    cudaFree(device_denominator_inv);
    cudaFree(device_eval_basis);
    cudaFree(device_merkle_root);
    cudaFree(device_layer_hashes);

    free(flattened_codeword);
    free(flattened_codeword_nxt);

    printf("Memory freed and commit_launch completed.\n");
}
