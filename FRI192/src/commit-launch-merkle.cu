
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "../include/commit-launch-merkle.cuh"
#include "../include/hash.cuh"
#include "../include/hash-host.cuh"
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

    uint64_t combined[2 * HASH_WORDS];
    uint64_t combined_codeword_hash[FIELD_WORDS + HASH_WORDS];
    uint8_t digest[HASH_SIZE];
    if (I < N / 2 && N == 131072) { //hardcoded for now
        memcpy(combined, &device_tree_layer[idx1], FIELD_WORDS * sizeof(uint64_t));  
        memcpy(combined + FIELD_WORDS, &device_tree_layer[idx2], FIELD_WORDS * sizeof(uint64_t));  
        SHA3(digest, (uint8_t *)combined, (FIELD_WORDS + FIELD_WORDS) * sizeof(uint64_t), 256);
        memcpy((uint8_t *)&device_layer_hashes[I * HASH_WORDS], digest, 32);
        memcpy(combined_codeword_hash, &device_codeword_nxt[idx3], FIELD_WORDS * sizeof(uint64_t));
        memcpy(combined_codeword_hash + FIELD_WORDS, digest, HASH_WORDS * sizeof(uint64_t));
        memcpy(&device_tree_layer_nxt[idx3], combined_codeword_hash, CONCAT_WORDS * sizeof(uint64_t));
    }
    __syncthreads();

    for (int layer = 1; layer < 12; layer++) { //12 is num_rounds
        if (I < (N / (2 * (1 << layer)))) {
            int codeword_idx = (1 << layer) * idx3;
            memcpy(combined, &device_layer_hashes[(I * 2) * HASH_WORDS], HASH_WORDS * sizeof(uint64_t));
            memcpy(combined + FIELD_WORDS, &device_tree_layer[codeword_idx], FIELD_WORDS * sizeof(uint64_t));

            SHA3(digest, (uint8_t *)combined,  FIELD_WORDS+HASH_WORDS * sizeof(uint64_t), 256);
            memcpy((uint8_t *)&device_layer_hashes[I * HASH_WORDS], digest, 32);
            memcpy(&device_tree_layer_nxt[I * CONCAT_WORDS], combined, CONCAT_WORDS * sizeof(uint64_t));
        }
        __syncthreads();
    }
}
__global__ void merkle_kernel(
    uint64_t *device_layer_hashes, 
    uint64_t *device_merkle_root, 
    uint64_t *device_tree_layer,
    uint64_t *device_tree_layer_nxt,
    int N
) {
    int I = blockIdx.x * blockDim.x + threadIdx.x;

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

        // // Store the combined hash values in `device_concatenated_tree` for this layer
        // int concat_index = I * (2 * HASH_WORDS);
        // memcpy(&device_concatenated_tree[concat_index], combined, 2 * HASH_WORDS * sizeof(uint64_t));
    }

    if (I == 0 && N == 2) {
        // Only copy to root when at the final layer
        memcpy(device_merkle_root, device_layer_hashes, HASH_WORDS * sizeof(uint64_t));
        printf("Final Merkle root: ");
        for (int j = 0; j < HASH_WORDS; j++) {
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
__global__ void final_layer_concat_and_hash_kernel(
    uint64_t *device_codeword_nxt,      // Device array of final layer codeword elements
    uint64_t *device_layer_hashes,      // Device array to store the updated hashes
    uint64_t *device_tree_layer, // Device array to store concatenated hash values
    uint64_t *device_tree_layer_nxt,
    int num_elements                    // Number of elements (N/2 == 32)
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < num_elements) {
        uint64_t combined[FIELD_WORDS + HASH_WORDS];
        uint8_t digest[HASH_SIZE];

        // Concatenate `codeword_nxt` with `layer_hashes`
        memcpy(combined, &device_codeword_nxt[i * FIELD_WORDS], FIELD_WORDS * sizeof(uint64_t));
        memcpy(combined + FIELD_WORDS, &device_layer_hashes[i * HASH_WORDS], HASH_WORDS * sizeof(uint64_t));

        // Hash the combined result
        SHA3(digest, (uint8_t *)combined, (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t), 256);

        // Copy the digest back into `device_layer_hashes`
        memcpy(&device_layer_hashes[i * HASH_WORDS], digest, HASH_SIZE);

        // Store the combined values in `device_concatenated_tree` for this layer
        memcpy(&device_tree_layer_nxt[i * CONCAT_WORDS], combined, (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
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

    if (N == 131072) {
        initialize_file("temp1.txt");
        initialize_file("temp2.txt");
        initialize_file("temp3.txt");
        initialize_file("temp4.txt");
        initialize_file("temp5.txt");
        initialize_file("alpha_offset.txt");
    }

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
    uint64_t *device_layer_hashes, *device_merkle_root, *device_tree_layer, *device_tree_layer_nxt;
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
    if(N == 131072){
    cudaMalloc((void**)&device_tree_layer, field_size);
    } else 
    {
        cudaMalloc((void**)&device_tree_layer, N * CONCAT_WORDS * sizeof(uint64_t));
    }
    cudaMalloc((void**)&device_tree_layer_nxt, N/2 * CONCAT_WORDS * sizeof(uint64_t));
    cudaMalloc((void**)&device_merkle_root, HASH_WORDS * sizeof(uint64_t));

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

//     // Variables to store memory usage
// size_t free_mem, total_mem;

// // Query the memory usage
// cudaMemGetInfo(&free_mem, &total_mem);

// // Calculate used memory
// size_t used_mem = total_mem - free_mem;

// // Print memory information
// printf("GPU Memory - Used: %zu MB, Free: %zu MB, Total: %zu MB\n", used_mem / (1024 * 1024), free_mem / (1024 * 1024), total_mem / (1024 * 1024));
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


    if (N == 64) { 
        int tree_idx = 2 * (int)log2(N);
        tree[tree_idx] = (uint64_t **)malloc((N/2) * sizeof(uint64_t *));
        for(int i = 0; i < N/2; i ++){
            tree[tree_idx][i] = (uint64_t *)malloc(CONCAT_WORDS * sizeof(uint64_t));
        }
        tree_idx++; 
        int next_N = N / 2;
        for(int layer = tree_idx; layer <= 17; layer++){
            tree[layer] = (uint64_t **)malloc((next_N / 2) * sizeof(uint64_t *));
            for(int i = 0; i < next_N / 2; i++){
                tree[layer][i] = (uint64_t *)malloc(2 * HASH_WORDS * sizeof(uint64_t));
                if (!tree[layer][i]) {
                    printf("Memory allocation failed at tree[%d][%d]\n", layer, i);
                    exit(1); 
                }
            }
            next_N = next_N / 2;
            printf("Allocated memory for tree[%d] with %d elements\n", layer, (next_N / 2) * HASH_WORDS);
        }
    
        int index = 13;
        int threads_per_block = 32;
        int num_blocks = (N / 2 + threads_per_block - 1) / threads_per_block;
    
        // Kernel for final concatenation and hash computation
        final_layer_concat_and_hash_kernel<<<num_blocks, threads_per_block>>>(
            device_codeword_nxt, 
            device_layer_hashes, 
            device_tree_layer,
            device_tree_layer_nxt,
            32
        );
        cudaDeviceSynchronize();
    
        // Copy the last flattened layer from device to host
        cudaMemcpy(flattened_tree_layer_nxt, device_tree_layer_nxt, (N / 2) * HASH_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    
        // Unflatten and assign to tree[12]
        for (int i = 0; i < N / 2; i++) {
            for (int j = 0; j < CONCAT_WORDS; j++) {
                tree[12][i][j] = flattened_tree_layer_nxt[i * CONCAT_WORDS + j];
            }
        }
    
        free(flattened_tree_layer_nxt);
        cudaFree(device_tree_layer_nxt);
    
        printf("Starting iterative Merkle kernel computation up the tree\n");
        N /= 2;
    
        // Allocate for iterative usage
        flattened_tree_layer_nxt = (uint64_t *)malloc((N / 2) * CONCAT_WORDS * sizeof(uint64_t));
        cudaMalloc((void **)&device_tree_layer_nxt, (N / 2) * CONCAT_WORDS * sizeof(uint64_t));
    
        while (N > 1) {
            int tpb = N / 2;
            int nb = (N + tpb - 1) / tpb;
    
            merkle_kernel<<<nb, tpb>>>(device_layer_hashes, device_merkle_root, device_tree_layer, device_tree_layer_nxt, N);
            cudaDeviceSynchronize();
    
            // Copy next layer from device to host
            cudaMemcpy(flattened_tree_layer_nxt, device_tree_layer_nxt, (N / 2) * HASH_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    
            for (int i = 0; i < N / 2; i++) {
                for (int j = 0; j < HASH_WORDS; j++) {
                    tree[index][i][j] = flattened_tree_layer_nxt[i * HASH_WORDS + j];
                }
            }
    
            if (N == 2) {
                cudaMemcpy(root, device_merkle_root, HASH_SIZE, cudaMemcpyDeviceToHost);
            }
    
            N /= 2;
            index++;
        }
    
        free(flattened_tree_layer_nxt);
        cudaFree(device_tree_layer_nxt);
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
