#ifndef FRI_CUDA_H
#define FRI_CUDA_H

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <cuda_runtime.h>  // CUDA Runtime
#include "../include/fiat-shamir.cuh"
#include "../include/merkle.cuh"
#include "../include/params.cuh"
#include "../include/field.cuh"
#include "../include/poly.cuh"
#include "../include/domain.cuh"
#include "../include/commit-launch-merkle.cuh"
#include "../include/poly-eval-launch.cuh"

int num_leaves = preon.fri_domains[5].size;
#define NUM_LEAVES num_leaves
#define MAX_PROOF_PATH_LENGTH ((int)log2(NUM_LEAVES)) //17
uint64_t *EVAL_BASIS_DICT[MAX_FRI_PARAMETERS][MAX_EVAL_BASIS_LEN];
uint64_t *OFFSET_DICT[MAX_FRI_PARAMETERS];
int basis_size = 0;

#define BASIS_SIZE 221
int push_count = 0;
int pull_count = 0;
uint64_t FRI_PARAMETERS[MAX_FRI_PARAMETERS][2][MAX_EVAL_BASIS_LEN];
uint64_t GF_2_192_BASIS[BASIS_SIZE * FIELD_WORDS];
const size_t field_words = 5;
extern uint64_t FRI_PARAMETERS[MAX_FRI_PARAMETERS][2][MAX_EVAL_BASIS_LEN];

typedef struct {
    int initial_domain_length;
    int expansion_factor;
    int num_colinearity_tests;
    int domain_length;
} Fri;

typedef struct {
    size_t *a_indices;
    size_t *b_indices;
    size_t *c_indices;
    int num_colinearity_tests;
} QueryIndices;


//extern functions
extern void hash_sha3_256(const uint64_t *data, size_t len, uint64_t *out);
extern void merkle_commit(uint64_t **leaf_hashes, size_t num_leaves, uint64_t *out);
extern void push_object(ProofStream *ps, void *obj);
extern unsigned char* prover_fiat_shamir(ProofStream *ps, size_t num_bytes);


__host__ __device__ void i_th_ele_in_span(uint64_t *result, uint64_t *basis, int len_basis, int i);
__host__ __device__ void field_addEqual(uint64_t *a, const uint64_t *b, size_t words);
__host__ __device__ void field_sub(uint64_t *out, const uint64_t *a, const uint64_t *b, size_t field_words);
__host__ __device__ void field_add(uint64_t *out, const uint64_t *a, const uint64_t *b, size_t field_words); 
__host__ __device__ void field_mul(uint64_t *out, const uint64_t *a, const uint64_t *b, size_t field_words);
__host__ __device__ void field_inv(uint64_t *out, const uint64_t *in, size_t field_words);
void extract_fri_parameters(const Domain *domains, int num_domains, uint64_t parameters[MAX_FRI_PARAMETERS][2][MAX_EVAL_BASIS_LEN]);
void init_fri_parameters(uint64_t parameters[MAX_FRI_PARAMETERS][2][MAX_EVAL_BASIS_LEN]);
void populate_eval_basis_and_offset_dict();
void cleanup_eval_basis_and_offset_dict();
void generate_elements(uint64_t *elements, const uint64_t *var1, const uint64_t *var2, size_t field_words);
Fri* init_fri(int initial_domain_length, int expansion_factor, int num_colinearity_tests);
int fri_num_rounds(Fri* fri);
__host__ __device__ void print_field(const char *label, const uint64_t *field, size_t field_words);

void byte_array_to_192bit_integer(const unsigned char *byte_array, size_t byte_array_len, uint64_t *output);
void field_sample(uint8_t *byte_array, size_t byte_array_len, uint64_t *eval_basis, size_t basis_len, uint64_t *result);
void int_to_bytes(int n, unsigned char *bytes, size_t size);
void sample_indices(uint8_t* seed, size_t seed_len, int size, int reduced_size, int number, size_t* indices, size_t* reduced_indices);
int fri_log_domain_length(Fri* fri);
void load_precomputed_inverses(const char *filename, uint64_t inverses[MAX_FRI_PARAMETERS]);
__host__ __device__ void print_codeword(uint64_t **codeword, size_t length, const size_t field_words);

void calculate_indices(uint64_t *c_indices, uint64_t *a_indices, uint64_t *b_indices, int num_colinearity_tests);
size_t* query(Fri *fri, uint64_t ***codewords, uint64_t **current_codeword, int current_codeword_len, uint64_t **next_codeword, size_t next_codeword_length, size_t *c_indices);
uint64_t* prove(Fri* fri, uint64_t **codeword);
int verify(Fri *fri, uint64_t **polynomial_values, int degree);

void print_eval_basis(const char *label, uint64_t *last_basis, size_t basis_len, size_t field_words);
void print_array(const char *label, const uint64_t *array, size_t size);

#endif /* FRI_CUDA_H */