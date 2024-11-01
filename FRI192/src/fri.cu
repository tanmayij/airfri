#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <cuda_runtime.h>
#include <math.h>
#include "../include/fiat-shamir.cuh"
#include "../include/fri.cuh"
#include "../include/merkle.cuh"
#include "../include/params.cuh"
#include "../include/field.cuh"
#include "../include/poly.cuh"
#include "../include/domain.cuh"
#include "../include/commit-launch-merkle.cuh"
#include "../include/poly-eval-launch.cuh"
#define FIELD_WORDS 4
ProofStream* global_proof_stream = NULL;
QueryIndices global_query_indices;

// typedef struct {
//     size_t size;
//     size_t basis_len;
//     const uint64_t *basis;
//     const uint64_t *shift;
// } Domain;
// Function to populate an eval basis array with non-zero basis elements for a given basis_len
uint64_t* populate_eval_basis(int basis_len) {
    int found = -1;
    // int max_basis_len = 80;
    for (int i = 0; i < sizeof(preon.fri_domains) / sizeof(preon.fri_domains[0]); i++) {
        if (preon.fri_domains[i].basis_len == basis_len) {
            found = i;
            break;
        }
    }

    if (found == -1) {
        fprintf(stderr, "Error: No FRI domain found with basis_len %d\n", basis_len);
        return NULL;
    }

    int count = 0;
    for (int i = 0; i < (preon.fri_domains[found].basis_len) * 4; i++) {
        if (preon.fri_domains[found].basis[i] != 0) {
            count++;
        }
    }
    uint64_t *eval_basis = (uint64_t *)malloc(count * sizeof(uint64_t));
    if (eval_basis == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for eval_basis\n");
        return NULL;
    }
    int index = 0;
    for (int i = 0; i < (preon.fri_domains[found].basis_len) * 4; i++) {
        if (preon.fri_domains[found].basis[i] != 0) {
            eval_basis[index++] = preon.fri_domains[found].basis[i];
        }
    }
    return eval_basis;
}

uint64_t get_offset_for_basis(int basis_len) {
    for (int i = 0; i < sizeof(preon.fri_domains) / sizeof(preon.fri_domains[0]); i++) {
        if (preon.fri_domains[i].basis_len == basis_len) {
            // Return the first non-zero element in the shift array for that domain
            for (int j = 0; j < 4; j++) {  //Assuming the shift array has 4 elements
                if (preon.fri_domains[i].shift[j] != 0) {
                    return preon.fri_domains[i].shift[j];
                }
            }
            fprintf(stderr, "Warning: No non-zero shift found for basis_len %d\n", basis_len);
            return 0;
        }
    }
    fprintf(stderr, "Error: No FRI domain found with basis_len %d\n", basis_len);
    return 0;
}

uint64_t* populate_all_bases(int *total_elements) {
    *total_elements = 0;
    //printf("tanjan %zu",sizeof(preon.fri_domains) / sizeof(preon.fri_domains[0]));
    // First pass to calculate the total number of non-zero elements dynamically
    for (size_t i = 0; i < sizeof(preon.fri_domains) / sizeof(preon.fri_domains[0]); i++) {
        for (int j = 0; j < (preon.fri_domains[i].basis_len) * 4; j++) {
            if (preon.fri_domains[i].basis[j] != 0ull) {
                (*total_elements)++;
            }
        }
    }
    uint64_t *all_bases = (uint64_t *)malloc(*total_elements * sizeof(uint64_t));
    if (all_bases == NULL) {
        fprintf(stderr, "Error: Memory allocation failed for all_bases\n");
        return NULL;
    }
    int basis_index = 0;
    for (int i = 0; i < sizeof(preon.fri_domains) / sizeof(preon.fri_domains[0]); i++) {
        for (int j = 0; j < (preon.fri_domains[i].basis_len) * 4; j++) {
            if (preon.fri_domains[i].basis[j] != 0ull) {
                all_bases[basis_index++] = preon.fri_domains[i].basis[j];
            }
        }
    }
    return all_bases; 
}

//the below function returns the ith element in the basis span
__host__ __device__ void i_th_ele_in_span(uint64_t *result, uint64_t *basis, int len_basis, int i)
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

//the below function generates elements by extending 2^64, based on field_words
void generate_elements(uint64_t *elements, const uint64_t *var1, const uint64_t *var2, const size_t field_words) {
    uint64_t /*replaced from uint32_t */ idx = 0;

    // Start with the elements 1, var1, var1^2, ..., var1^63
    uint64_t temp[field_words];
    field_one(&elements[idx * field_words], field_words);
    idx++;

    for (int i = 1; i <= 63; i++) {
        field_pow(temp, var1, i, field_words);
        memcpy(&elements[idx * field_words], temp, field_words * sizeof(uint64_t));
        idx++;
    }

    // Now add var2, var1*var2, var1^2*var2, ..., var1^63*var2
    for (int i = 0; i <= 63; i++) {
        field_pow(temp, var1, i, field_words);
        field_mul(temp, temp, var2, field_words);
        memcpy(&elements[idx * field_words], temp, field_words * sizeof(uint64_t));
        idx++;
    }

    //now, add var2^2, var1*var2^2, var1^2*var2^2, ..., var1^63*var2^2
    uint64_t var2_squared[field_words];
    field_mul(var2_squared, var2, var2, field_words);

    for (int i = 0; i <= 63; i++) {
        field_pow(temp, var1, i, field_words);
        field_mul(temp, temp, var2_squared, field_words);
        memcpy(&elements[idx * field_words], temp, field_words * sizeof(uint64_t));
        idx++;
    }
    if(field_words == 4){
    //finally add var2^3, var1*var2^3, var1^2*var2^3, ..., var1^63*var2^3
    uint64_t var2_cube[field_words];
    field_mul(var2_cube, var2_squared, var2, field_words);

    for (int i = 0; i <= 63; i++) {
        field_pow(temp, var1, i, field_words);
        field_mul(temp, temp, var2_cube, field_words);
        memcpy(&elements[idx * field_words], temp, field_words * sizeof(uint64_t));
        idx++;
    }
    }
    if(field_words == 5){
        //finally add var2^3, var1*var2^3, var1^2*var2^3, ..., var1^63*var2^3
        uint64_t var2_four[field_words];
        uint64_t var2_cube[field_words];
        field_mul(var2_cube, var2_squared, var2, field_words);
        field_mul(var2_four, var2_cube, var2, field_words);
    
        for (int i = 0; i <= 63; i++) {
            field_pow(temp, var1, i, field_words);
            field_mul(temp, temp, var2_four, field_words);
            memcpy(&elements[idx * field_words], temp, field_words * sizeof(uint64_t));
            idx++;
        }
    }
}

//FRI initialization
Fri* init_fri(int initial_domain_length, int expansion_factor, int num_colinearity_tests) {
    Fri* fri = (Fri*)malloc(sizeof(Fri));
    fri->initial_domain_length = initial_domain_length;
    fri->expansion_factor = expansion_factor;
    fri->num_colinearity_tests = num_colinearity_tests;
    return fri;
}

//the below function computes the number of rounds in the protocol. it should be 12 consistently. 
int fri_num_rounds(Fri* fri) {
    int codeword_length = fri->initial_domain_length;
    printf("Initial codeword length: %d\n", codeword_length);

    int num_rounds = 0;
    int max_rounds = 12;  //set the maximum number of rounds to 12

    //loop until we reach 12 rounds or the codeword length becomes too small
    while (num_rounds < max_rounds && codeword_length > fri->expansion_factor && fri->num_colinearity_tests < codeword_length) {
        codeword_length /= 2;
        num_rounds += 1;
    }

    printf("Final codeword length after %d rounds: %d\n", num_rounds, codeword_length);
    return num_rounds;
}

//the below function prints field elements for debugging
__host__ __device__ void print_field(const char *label, const uint64_t *field, size_t field_words) {
    printf("%s: ", label);
    for (size_t i = 0; i < field_words; i++) {
        printf("%016lx ", field[i]);
    }
    printf("\n");
}

//function to copy field elements
void field_copy(uint64_t *dest, const uint64_t *src, size_t field_words) {
    memcpy(dest, src, field_words * sizeof(uint64_t));
}

//the below function samples elements from extension field based on the eval_basis. this specific function below works for 256-bit elements
void field_sample(uint8_t *byte_array, size_t byte_array_len, uint64_t *eval_basis, size_t basis_len, uint64_t *result) {
    uint64_t acc[4] = {0};  // 256 bits (4 * 64 bits)

    printf("Byte array: ");
    for (size_t i = 0; i < byte_array_len; ++i) {
        printf("%02x", byte_array[i]);
    }
    printf("\n");

    for (uint64_t i = 0; i < byte_array_len; ++i) {
        for (uint64_t j = 3; j > 0; --j) {
            acc[j] = (acc[j] << 8) | ((acc[j-1] >> 56) & 0xFF);
        }
        acc[0] = (acc[0] << 8) | byte_array[i];
    }

    printf("Accumulator after conversion: ");
    for (uint64_t i = 0; i < 4; i++) {
        printf("%016lx ", acc[i]);
    }
    printf("\n");

    //mask to ensure acc is at most 256 bits (though acc is already 256 bits, this is just to be sure)
    acc[0] &= 0xFFFFFFFFFFFFFFFF;
    acc[1] &= 0xFFFFFFFFFFFFFFFF;
    acc[2] &= 0xFFFFFFFFFFFFFFFF;
    acc[3] &= 0xFFFFFFFFFFFFFFFF;

    printf("Accumulator after masking: ");
    for (uint64_t i = 0; i < 4; i++) {
        printf("%016lx ", acc[i]);
    }
    printf("\n");

    uint64_t combined_index = acc[0] ^ acc[1] ^ acc[2] ^ acc[3];
    combined_index = combined_index % basis_len;  
    printf("Combined index: %016lx\n", combined_index);
    i_th_ele_in_span(result, eval_basis, basis_len, combined_index);

    //print the result for debugging
    printf("Sampled field element: ");
    for (uint64_t i = 0; i < FIELD_WORDS; i++) {
        printf("%016lx ", result[i]);
    }
    printf("\n");
}

//the below function sampled indices for use in query function
void sample_indices(uint8_t* seed, size_t seed_len, int size, int reduced_size, int number, size_t* indices, size_t* reduced_indices) {
    if (number > reduced_size) {
        fprintf(stderr, "cannot sample more indices than available in last codeword; requested: %d, available: %d\n", number, reduced_size);
        exit(EXIT_FAILURE);
    }
    if (number > 2 * reduced_size) {
        fprintf(stderr, "not enough entropy in indices wrt last codeword\n");
        exit(EXIT_FAILURE);
    } 
    printf("Seed in sample_indices: ");
    for (size_t i = 0; i < seed_len; i++) {
        printf("%02x", seed[i]);
    }
    printf("\n");
    printf("size: %d\n", size);
    printf("reduced_size: %d\n", reduced_size);
    int counter = 0;
    int indices_count = 0;
    int reduced_indices_count = 0;
    while (indices_count < number) { //indices_count < 14
        unsigned char counter_bytes[15];
        int_to_bytes(counter, counter_bytes, sizeof(counter_bytes));

        unsigned char *concatenated = (unsigned char *)malloc(seed_len + sizeof(counter_bytes));
        memcpy(concatenated, seed, seed_len);
        memcpy(concatenated + seed_len, counter_bytes, sizeof(counter_bytes));

        uint64_t digest[HASH_WORDS];
        hash_sha3_256((uint64_t *)concatenated, seed_len + sizeof(counter_bytes), digest);

        size_t index = 0;
        for (int i = 0; i < HASH_SIZE; i++) {
            index = (index << 8) | ((unsigned char *)digest)[i];
        }
        // printf("Index before modulo: %zu, size: %d\n", index, size);
        if (size <= 0) {
            printf("Error: size is invalid (%d), skipping modulo.\n", size);
            exit(EXIT_FAILURE);
        }
        index = index % size;
        //printf("Index after modulo: %lu\n", index);

        size_t reduced_index = index >> (int)(log2(size / reduced_size));

        counter += 1;  
        int found = 0;
        for (size_t i = 0; i < reduced_indices_count; i++) {
            if (reduced_indices[i] == reduced_index) {
                found = 1;
                break;
            }
        }

        if (!found) {
            indices[indices_count++] = index;
            reduced_indices[reduced_indices_count++] = reduced_index;
        }

        free(concatenated);  // Free allocated memory
    }

    // Sort indices
    for (size_t i = 0; i < indices_count - 1; i++) {
        for (size_t j = i + 1; j < indices_count; j++) {
            if (indices[i] > indices[j]) {
                size_t temp = indices[i];
                indices[i] = indices[j];
                indices[j] = temp;
            }
        }
    }

    // Sort reduced_indices
    for (size_t i = 0; i < reduced_indices_count - 1; i++) {
        for (size_t j = i + 1; j < reduced_indices_count; j++) {
            if (reduced_indices[i] > reduced_indices[j]) {
                size_t temp = reduced_indices[i];
                reduced_indices[i] = reduced_indices[j];
                reduced_indices[j] = temp;
            }
        }
    }

    // Print indices
    printf("Indices: ");
    for (size_t i = 0; i < indices_count; i++) {
        printf("%zu ", indices[i]);
    }
    printf("\n");

    // Print freduced_indices
    printf("Reduced Indices: ");
    for (size_t i = 0; i < reduced_indices_count; i++) {
        printf("%zu ", reduced_indices[i]);
    }
    printf("\n");
}

int fri_log_domain_length(Fri* fri) {
    return (int)log2(fri->initial_domain_length);
}
//below function reads precomputed inverses from the file and puts it in the array that the function call is assigned to  
void load_precomputed_inverses(const char *filename, uint64_t inverses[MAX_FRI_PARAMETERS]) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    char buffer[1024];  // To hold the lines from the file
    for (int i = 0; i < MAX_FRI_PARAMETERS; i++) {
        // Read lines until you find a line starting with "FRI Domain"
        while (fgets(buffer, sizeof(buffer), file)) {
            if (strstr(buffer, "FRI Domain") != NULL) {
                break;  // We found the correct domain, so exit this loop
            }
        }

        // Now read the first inverse element
        if (fscanf(file, "%016lx", &inverses[i]) != 1) {
            fprintf(stderr, "Error reading inverse for FRI Domain %d\n", i);
            fclose(file);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
}

//prints the whole codeword- for debugging 
__host__ __device__ void print_codeword(uint64_t **codeword, size_t length, const size_t field_words) {
    for (size_t i = 0; i < length; i++) {
        printf("codeword[%zu /*replaced from uint32_t */]: ", i);
        for (size_t j = 0; j < field_words; j++) {
            printf("%016lx ", codeword[i][j]);
        }
        printf("\n");
    }
}

//struct handles indices because there should not be an inconsistency in indices between query and verify. 
typedef struct {
    size_t **round_indices;  //array of pointers to indices for each round
    int num_rounds;           //number of FRI rounds
    int num_tests;            
} IndicesTracker;

IndicesTracker global_indices_tracker;
void initialize_indices_tracker(Fri* fri) {
    int num_rounds = fri_num_rounds(fri);
    global_indices_tracker.num_rounds = num_rounds;
    global_indices_tracker.num_tests = fri->num_colinearity_tests;

    // Allocate memory to store indices for each round
    global_indices_tracker.round_indices = (size_t **)malloc(num_rounds * sizeof(size_t *));
    for (int i = 0; i < num_rounds; i++) {
        global_indices_tracker.round_indices[i] = (size_t *)malloc(global_indices_tracker.num_tests * sizeof(size_t));
    }
}

void cleanup_indices_tracker() {
    for (int i = 0; i < global_indices_tracker.num_rounds; i++) {
        free(global_indices_tracker.round_indices[i]);
    }
    free(global_indices_tracker.round_indices);
}

/* ======== COMMIT PHASE ======== */

uint64_t ***commit_host(Fri *fri, uint64_t **codeword, int codeword_len)
{
    const int max_fri_domains = sizeof(preon.fri_domains) / sizeof(preon.fri_domains[0]); //should be 16 for 256-bit params 
    int basis_len = fri_log_domain_length(fri); //basis length of initial codeword 
    // uint64_t *eval_basis = populate_eval_basis(basis_len);
    // if (eval_basis != NULL) {
    //     printf("Eval basis in commit phase for basis_len %d:\n", basis_len);
    //     for (int i = 0; i < basis_len; i++) {
    //         printf("%016lx ", eval_basis[i]);
    //     }
    //     printf("\n");
    // }
    // assert(eval_basis != NULL && "eval_basis are not correct");
    // uint64_t offset = get_offset_for_basis(basis_len);

    // if (offset != 0) {
    //     printf("Offset in commit phase for basis_len %d is: %016lx\n", basis_len, offset);
    // } else {
    //     printf("Offset in commit phase not found or is zero.\n");
    // } 
    uint64_t ***codewords = (uint64_t ***)malloc((fri_num_rounds(fri) + 1 /* 12 + 1 */) * sizeof(uint64_t **)); //because we want 12 codewords.
    uint64_t precomputed_inverses[max_fri_domains];
    load_precomputed_inverses("fri_inverses_256.txt", precomputed_inverses);
    
    uint64_t var1[FIELD_WORDS] = {0x2}; // Represents var1 for generating elements in the field 
    uint64_t var2[FIELD_WORDS] = {0x3}; // Represents var2 for generating elements in the field 
    uint64_t elements[256 * FIELD_WORDS];
    generate_elements(elements, var1, var2, field_words);

    //below code samples alpha
    uint64_t *alpha = (uint64_t *) malloc (FIELD_WORDS * sizeof(uint64_t));
    uint64_t *sampled_alpha = (uint64_t *) malloc (FIELD_WORDS * sizeof(uint64_t));
    size_t N = codeword_len;
    printf("num rounds in commit are %d\n", fri_num_rounds(fri));
    size_t num_butes_alpha = 24;
    uint8_t *alpha_bytes = prover_fiat_shamir(global_proof_stream, num_butes_alpha); //initialise alpha_bytes
    field_sample(alpha_bytes, 24, elements, 256, sampled_alpha);
    printf("Sampled Alpha in commit: ");
    for (size_t i = 0; i < FIELD_WORDS; i++) {
        printf("%016lx ", sampled_alpha[i]);
    }
    printf("\n");
    uint64_t *sampled_alpha_for_proofstream = (uint64_t *) malloc (FIELD_WORDS * sizeof(uint64_t));
    memcpy(sampled_alpha_for_proofstream, sampled_alpha, FIELD_WORDS * sizeof(uint64_t));

    uint64_t *root = (uint64_t *)malloc(4 * sizeof(uint64_t)); //32 bytes for merkle tree hash
    for(int r = 0; r < fri_num_rounds(fri); r++){
        uint64_t *eval_basis = populate_eval_basis(basis_len - r);
        if (eval_basis != NULL) {
            printf("Eval basis for basis_len %d:\n", basis_len);
            for (int i = 0; i < basis_len; i++) {
                printf("%016lx ", eval_basis[i]);
            }
            printf("\n");
        }
        printf("round: %d\n", r);
        uint64_t offset = get_offset_for_basis(basis_len - r);
        if (offset != 0) {
            printf("Offset for basis_len %d is: %016lx\n", basis_len, offset);
        } else {
            printf("Offset not found or is zero.\n");
        }
        assert(eval_basis!=NULL && "eval_basis in commit are not correct");

        //hash alpha repeatedly for every round
        hash_sha3_256((uint64_t *) sampled_alpha, 24, (uint64_t *) alpha);
        memcpy(sampled_alpha, alpha, 24);
        printf("Alpha that is hashed in commit: ");
        for (size_t i = 0; i < FIELD_WORDS; i++) {
            printf("%016lx ", alpha[i]);
        }
        printf("\n");

        codewords[r] = codeword;
        int exponent = preon.fri_domains[0].basis_len - fri_log_domain_length(fri);
        printf("preon.fri_domains[0].basis_len: %d\n", preon.fri_domains[0].basis_len);
        uint64_t denominator_inv = 0;
        denominator_inv = precomputed_inverses[exponent + r]; 
        printf("denominator inverse %016lx\n", denominator_inv);
        uint64_t **codeword_nxt = (uint64_t **)malloc((N/2)*sizeof(uint64_t *)); 

        commit_launch(codeword, codeword_nxt, alpha, &offset, denominator_inv, eval_basis, N, root);
        if(N > 32) {
        N = N / 2;
        codeword = codeword_nxt;
        codeword_len = codeword_len /2;
        }
    }
    push_object(global_proof_stream, root);
    push_count++;
    push_object(global_proof_stream, sampled_alpha_for_proofstream);
    push_count++;
    printf("push counted +1 in commit for root and sampled alpha %d", push_count);
    printf("\n");
    printf("Root of FRI Merkle tree:\n");
    for (size_t i = 0; i < 4; i++) {
        printf("%016lx ", root[i]);
    }
    printf("\n");
    printf("Push successful\n");
    //push the last codeword and store it
    for (int i =0; i < N; i ++) {
    push_object(global_proof_stream, codeword[i]);
    push_count++;
    }
    codewords[fri_num_rounds(fri)] = codeword;

    //debugging information
    print_codeword(codewords[fri_num_rounds(fri)], N, field_words);
    printf("Last codeword length N:%zu", N);

    //return the last codeword
    return codewords;
}

/* ======== QUERY PHASE ======== */

void calculate_indices(size_t  *c_indices, size_t  *a_indices, size_t  *b_indices, int num_colinearity_tests) {

    for (int i = 0; i < num_colinearity_tests; i++) {
        a_indices[i] = (2 * c_indices[i]);
        b_indices[i] = (2 * c_indices[i] + 1);
    }
}

size_t* query(Fri *fri, uint64_t ***codewords, uint64_t **current_codeword, size_t current_codeword_len, uint64_t **next_codeword, size_t next_codeword_len, size_t *c_indices) 
{
    int num_tests = fri->num_colinearity_tests;
    global_query_indices.a_indices = (size_t *)malloc(num_tests * sizeof(size_t));
    global_query_indices.b_indices = (size_t *)malloc(num_tests * sizeof(size_t));
    global_query_indices.c_indices = (size_t *)malloc(num_tests * sizeof(size_t));
    global_query_indices.num_colinearity_tests = num_tests;
    memcpy(global_query_indices.c_indices, c_indices, global_query_indices.num_colinearity_tests * sizeof(size_t));
    // Calculate a and b indices
    calculate_indices(global_query_indices.c_indices, 
    global_query_indices.a_indices, 
    global_query_indices.b_indices, global_query_indices.num_colinearity_tests);

    //the below loop is going to run for each round, for the total number of tests. 
    for(int s = 0; s < num_tests; s++){
        push_object(global_proof_stream, current_codeword[global_query_indices.a_indices[s]]);
        print_field("current_codeword[a_indices[s]]", current_codeword[global_query_indices.a_indices[s]], field_words);
        push_count++;
        push_object(global_proof_stream, current_codeword[global_query_indices.b_indices[s]]);
        print_field("current_codeword[b_indices[s]]", current_codeword[global_query_indices.b_indices[s]], field_words);
        push_count++;
        push_object(global_proof_stream, next_codeword[global_query_indices.c_indices[s]]);
        print_field("next_codeword[c_indices[s]]", next_codeword[global_query_indices.c_indices[s]], field_words);
        push_count++;

        //below logic is to compute auth_paths for each of the codewords and their indices.

        //for a_index
        size_t proof_len_a = 0;
        uint64_t **auth_path_a = (uint64_t **)malloc(MAX_PROOF_PATH_LENGTH * sizeof(uint64_t *));
        for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
        auth_path_a[i] = (uint64_t *)malloc(HASH_SIZE);  // Allocate space for each hash
        }
        merkle_open(codewords, current_codeword_len, global_query_indices.a_indices[s], auth_path_a, &proof_len_a, field_words);
        push_object(global_proof_stream, &proof_len_a);
        push_count++;
        for (size_t i = 0; i < proof_len_a; i++) {
        //auth_path_a[%zu]\n", i);
        push_object(global_proof_stream, auth_path_a[i]); // Push each hash
        push_count++;
        }
        printf("proof len a: %zu \n", proof_len_a);
        printf("push counted +1 in query for a_indices %d\n", push_count);
        for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
            free(auth_path_a[i]);
        }
        free(auth_path_a);

        //for b_index
        size_t proof_len_b = 0;
        uint64_t **auth_path_b = (uint64_t **)malloc(MAX_PROOF_PATH_LENGTH * sizeof(uint64_t *));
        for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
        auth_path_b[i] = (uint64_t *)malloc(HASH_SIZE);  //allocate space for each hash
        }
        merkle_open(codewords, current_codeword_len, global_query_indices.b_indices[s], auth_path_b, &proof_len_b, field_words);
        push_object(global_proof_stream, &proof_len_b);
        push_count++;
        for (size_t i = 0; i < proof_len_b; i++) {
        //printf("Pushing auth_path_b[%zu]\n", i);
        push_object(global_proof_stream, auth_path_b[i]); //push each hash
        push_count++;
        }
        printf("proof len b: %zu \n", proof_len_b);
        printf("push counted +1 in query for b_indices %d\n", push_count);
        for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
            free(auth_path_b[i]);
        }
        free(auth_path_b);

        //for c_index
        size_t proof_len_c = 0;
        uint64_t **auth_path_c = (uint64_t **)malloc(MAX_PROOF_PATH_LENGTH * sizeof(uint64_t *));
        for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
        auth_path_c[i] = (uint64_t *)malloc(HASH_SIZE);  //allocate space for each hash
        }
        merkle_open(codewords, next_codeword_len, global_query_indices.c_indices[s], auth_path_c, &proof_len_c, field_words);
        push_object(global_proof_stream, &proof_len_c);
        push_count++;
        for (size_t i = 0; i < proof_len_c; i++) {
        //printf("Pushing auth_path_c[%zu]\n", i);
        push_object(global_proof_stream, auth_path_c[i]); //push each hash
        push_count++;
        }
        printf("proof len c: %zu \n", proof_len_c);
        for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
            free(auth_path_c[i]);
        }
        free(auth_path_c);
    }
    printf("Indices in query:\n");
    for (int i = 0; i < fri->num_colinearity_tests; i++) {
    printf("a: %zu, b: %zu, c: %zu\n", global_query_indices.a_indices[i], global_query_indices.b_indices[i], global_query_indices.c_indices[i]);
    }
    return c_indices;
}

size_t* prove(Fri* fri, uint64_t **codeword) {
    initialize_indices_tracker(fri);

    size_t codeword_length = 0;
    while (codeword[codeword_length] != NULL) codeword_length++;

    global_proof_stream = create_proof_stream();
    uint64_t ***codewords = commit_host(fri, codeword, codeword_length);

    size_t num_bytes = 32; 
    uint8_t *seed = prover_fiat_shamir(global_proof_stream, num_bytes); 
    //Initialize space for top-level and reduced indices
    size_t *top_level_indices = (size_t *)malloc(fri->num_colinearity_tests * sizeof(size_t));
    size_t *reduced_indices = (size_t *)malloc(fri->num_colinearity_tests * sizeof(size_t));

    //sample indices for the first round
    sample_indices(seed, 32, codeword_length / 2, (fri->initial_domain_length >> (fri_num_rounds(fri))), fri->num_colinearity_tests, top_level_indices, reduced_indices);
    size_t *indices = (size_t *)malloc(fri->num_colinearity_tests * sizeof(size_t));
    memcpy(indices, top_level_indices, fri->num_colinearity_tests * sizeof(size_t));

    for (size_t i = 0; i < fri_num_rounds(fri) - 1; i++) { 
        query(fri, codewords, codewords[i], codeword_length, codewords[i + 1], codeword_length / 2, indices);
        memcpy(global_indices_tracker.round_indices[i], indices, fri->num_colinearity_tests * sizeof(size_t));

        for (size_t j = 0; j < fri->num_colinearity_tests; j++) {
            indices[j] /= 2;
        }
        codeword_length /= 2;
    }

    free(seed);
    free(top_level_indices);
    free(reduced_indices);
    free(indices);

    return top_level_indices;
}

void print_element(uint64_t *element, size_t num_words) {
    for (size_t i = 0; i < num_words; i++) {
        printf("%016lx ", element[i]);
    }
    printf("\n");
}


int verify(Fri *fri, uint64_t **polynomial_values, int degree) {
    uint64_t one[FIELD_WORDS];
    field_one(one, FIELD_WORDS);
    const int max_fri_domains = sizeof(preon.fri_domains) / sizeof(preon.fri_domains[0]); //should be 16 for 256-bit params 
    int basis_len = fri_log_domain_length(fri); //basis length of initial codeword 
    uint64_t *eval_basis = populate_eval_basis(basis_len);
    if (eval_basis != NULL) {
        printf("Eval basis in commit phase for basis_len %d:\n", basis_len);
        for (int i = 0; i < basis_len; i++) {
            printf("%016lx ", eval_basis[i]);
        }
        printf("\n");
    }
    assert(eval_basis != NULL && "eval_basis are not correct");
    uint64_t offset = get_offset_for_basis(basis_len);

    if (offset != 0) {
        printf("Offset in commit phase for basis_len %d is: %016lx\n", basis_len, offset);
    } else {
        printf("Offset in commit phase not found or is zero.\n");
    } 
    assert(eval_basis != NULL && "eval_basis are not correct!");
    uint64_t *root_verify = (uint64_t *)malloc(HASH_SIZE * sizeof(uint64_t)); 
    size_t last_codeword_length = (int)pow(2, fri_log_domain_length(fri)) >> fri_num_rounds(fri); //should be 32 for FRI192

    //get the elements here too
    uint64_t var1[FIELD_WORDS] = {0x2}; // Represents var1
    uint64_t var2[FIELD_WORDS] = {0x3}; // Represents var2
    uint64_t elements[256*FIELD_WORDS];

    generate_elements(elements, var1, var2, field_words);
    uint64_t *alpha = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
    uint64_t *sampled_alpha = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
    uint8_t *alpha_bytes = verifier_fiat_shamir(global_proof_stream, 24);

    root_verify = (uint64_t *)pull_object(global_proof_stream);
    pull_count++;
    sampled_alpha = (uint64_t *)pull_object(global_proof_stream);
    pull_count++;

    printf("printing FRI Merkle root");
    for(int j = 0; j< FIELD_WORDS; j++){
        printf("%016lx ", root_verify[j]);
    }
    printf("\n");

    //get the last codeword elements 
    printf("last codeword length: %zu\n", last_codeword_length);
    uint64_t **last_codeword_arr = (uint64_t **)malloc(last_codeword_length * sizeof(uint64_t));
    for (int i = 0; i < last_codeword_length; i ++){
        *last_codeword_arr = (uint64_t *)malloc(field_words * sizeof(uint64_t *));
    }   
    //extract last codeword
    for(int i=0; i<last_codeword_length; i++){
        last_codeword_arr[i] = (uint64_t *)pull_object(global_proof_stream);
        pull_count++;
    }
    printf("\n");

    uint64_t *copied_last_codeword = (uint64_t *)malloc(HASH_SIZE);
    memcpy(copied_last_codeword, last_codeword_arr, HASH_SIZE);

    //degree fetch
    degree = (last_codeword_length * fri_num_rounds(fri) / fri->expansion_factor) - 1; //if last codeword length = 32, this should be 0
    //last basis and offset
    uint64_t *last_basis = (uint64_t *)malloc((int)log2(last_codeword_length) * sizeof(uint64_t));
    last_basis = populate_eval_basis((int)log2(last_codeword_length));
    uint64_t last_offset = get_offset_for_basis((int)log2(last_codeword_length));

    Domain last_domain;
    int r = fri_num_rounds(fri);

    initialize_domain(&last_domain, fri->initial_domain_length >> r, last_basis, &last_offset);
    uint64_t **last_domain_elements = (uint64_t **)malloc(last_domain.size * sizeof(uint64_t *));
    for (int i = 0; i < last_domain.size; i++) {
        last_domain_elements[i] = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t)); //this is the interpolant
        i_th_ele_in_span(last_domain_elements[i], last_basis, last_domain.basis_len, i);
        field_add(last_domain_elements[i], last_domain_elements[i], &last_offset, 1);
    }

    uint64_t poly_coeffs_[field_words * last_domain.size];
    interpolate_domain_single(poly_coeffs_, last_domain_elements, copied_last_codeword, last_domain.size, FIELD_WORDS);
    uint64_t eval_result[FIELD_WORDS];
    poly_eval_over_domain(eval_result, poly_coeffs_, last_domain.size, &last_domain, field_words, preon.field_bytesize); //this function calls fft
    if (memcmp(copied_last_codeword, eval_result, FIELD_WORDS * sizeof(uint64_t)) == 0) {
        printf("Re-evaluated codeword does not match original!\n");
        free(root_verify);
        for (int j = 0; j < last_domain.size; j++) {
            free(last_domain_elements[j]);
        }
        free(last_domain_elements);
        return 0;
    }
    if (poly_deg(poly_coeffs_, last_domain.size, FIELD_WORDS) <= degree) {
        printf("Low degree test passed\n");
        printf("Observed degree: %d\n", poly_deg(poly_coeffs_, last_domain.size, FIELD_WORDS) );
        printf("Should be no more than: %d\n", degree);
    } else {
        printf("Last codeword does not correspond to polynomial of low enough degree\n");
        printf("Observed degree: %d\n", poly_deg(poly_coeffs_, last_domain.size, FIELD_WORDS));
        printf("But should be no more than: %d\n", degree);
        return 0;
    } 
    for (r = 0; r < fri_num_rounds(fri); r++) {
        const int max_fri_domains = sizeof(preon.fri_domains) / sizeof(preon.fri_domains[0]);
        int basis_len = fri_log_domain_length(fri) >> (r); 
        eval_basis = populate_eval_basis(basis_len);
        offset = get_offset_for_basis(basis_len);
        assert(eval_basis != NULL && "eval_basis are not correct!");

        //hash alpha
        hash_sha3_256((uint64_t *) sampled_alpha, 24, (uint64_t *) alpha); //hash alpha
        memcpy(sampled_alpha, alpha, 24); //repeat hash
        printf("Alpha that is hashed in verify: ");
        for (uint64_t  i = 0; i < field_words; i++) {
            printf("%016lx ", alpha[i]);
        }
        printf("\n");

        //get the indices
        size_t *c_indices = global_indices_tracker.round_indices[r];
        calculate_indices(c_indices, 
            global_query_indices.a_indices, 
            global_query_indices.b_indices, 
            fri->num_colinearity_tests);
        
        printf("Indices in verify:\n");
        for (int i = 0; i < fri->num_colinearity_tests; i++) {
            printf("a: %zu, b: %zu, c: %zu\n", global_query_indices.a_indices[i], global_query_indices.b_indices[i], c_indices[i]);
        }

        //arrays store information about codeword elements for each colinearity test
        uint64_t **aa = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **bb = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **cc = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **ay = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **by = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **cy = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));

        uint64_t precomputed_inverses[max_fri_domains];
        load_precomputed_inverses("fri_inverses_256.txt", precomputed_inverses);
        int exponent = preon.fri_domains[0].basis_len - fri_log_domain_length(fri); 
        uint64_t denominator_inv[field_words] = {0};
    
        denominator_inv[0] = precomputed_inverses[exponent + r]; 
        printf("denominator inverse %016lx\n", denominator_inv[0]);
        denominator_inv[1] = {0};
        denominator_inv[2] = {0};
        denominator_inv[3] = {0};
        //colinearity check. to change: ideally, we must check if merkle root verifies and then perform check. but for now this order is okay.
        for(int s = 0; s < fri -> num_colinearity_tests; s++){
            ay[s] = (uint64_t *)pull_object(global_proof_stream);
            pull_count++;
            print_field("ay", ay[s], field_words);
            aa[s] = ay[s];
            by[s] = (uint64_t *)pull_object(global_proof_stream);
            pull_count++;
            print_field("by", by[s], field_words);
            bb[s] = by[s];
            cy[s] = (uint64_t *)pull_object(global_proof_stream);
            pull_count++;
            print_field("cy", cy[s], field_words);
            cc[s] = cy[s];
         
            if (r == 0) {
                polynomial_values[global_query_indices.a_indices[s]] = aa[s];
                polynomial_values[global_query_indices.b_indices[s]] = bb[s];
            }
            // Printing the polynomial values that are not zero
            printf("Polynomial values:\n");
            int k = 0; // Define k outside the loop
            for (int i = 0; i < fri->num_colinearity_tests; i++) {
                //check if the polynomial value is zero, skip if true
                if (is_zero(polynomial_values[i], 1)) {
                    continue;
                }
            
                // Print the current non-zero polynomial value
                for (int j = 0; j < FIELD_WORDS; j++) {  //assuming FIELD_WORDS is the length of each polynomial value
                    printf("%016lx ", polynomial_values[i][j]);
                }
                printf("\n");
                for (int j = 0; j < FIELD_WORDS; j++) {
                    polynomial_values[k][j] = polynomial_values[i][j];
                }
                k++;
            }
            uint64_t ax[field_words] = {0};
            uint64_t ax_temp[field_words] = {0};
            uint64_t bx[field_words] = {0};
            uint64_t bx_temp[field_words] = {0};
            uint64_t alpha_offset[field_words] = {0};
            uint64_t cx[field_words] = {0};
            uint64_t diff1[field_words] = {0};
            uint64_t diff2[field_words] = {0};
            uint64_t temp1[field_words] = {0};
            uint64_t temp2[field_words] = {0};
            uint64_t temp3[field_words] = {0};
            uint64_t temp4[field_words] = {0};
            uint64_t temp5[field_words] = {0};

            i_th_ele_in_span(ax_temp, eval_basis, basis_len - r, global_query_indices.a_indices[s]);
            field_add(ax, ax_temp, &offset, field_words);
            print_field("ax", ax, field_words);
            ax[1] = {0};
            ax[2] = {0};
            ax[3] = {0};
            i_th_ele_in_span(bx_temp, eval_basis, basis_len - r, global_query_indices.b_indices[s]);
            field_add(bx, &offset, bx_temp, field_words);
            print_field("bx", bx, field_words);
            bx[1] = {0};
            bx[2] = {0};
            bx[3] = {0};
            memcpy(cx, alpha, FIELD_WORDS * sizeof(uint64_t));
            field_sub(diff1, ay[s], by[s], field_words);
            print_field("diff1", diff1, field_words);
            field_sub(diff2, ax, bx, field_words);
            field_inv(temp1, diff2, field_words);
            memcpy(alpha_offset, ax, sizeof(uint64_t));
            field_sub(temp2, cx, ax, field_words);
            print_field("cx", cx, field_words);
            print_field("temp2", temp2, field_words);
            field_mul(temp3, temp2, denominator_inv, field_words);
            print_field("temp3", temp3, field_words);
            field_mul(temp4, temp3, diff1, field_words);
            print_field("temp4", temp4, field_words);
            field_add(temp5, temp4, ay[s], field_words);
            print_field("temp5", temp5, field_words);
            if (memcmp(cy[s], temp5, field_words * sizeof(uint64_t))!= 0) {
                printf("Colinearity check failed\n");
                printf("%d\n", memcmp(cy[s], temp5, field_words * sizeof(uint64_t)));
                free(aa);
                free(bb);
                free(cc);
                //free(root_verify);
                //free(alpha);
                for (int j = 0; j < last_domain.size; j++) {
                    free(last_domain_elements[j]);
                }
                free(last_domain_elements);
                return 0;
            }

            //check auth_path_a
            uint64_t **auth_path_a = (uint64_t **)malloc(MAX_PROOF_PATH_LENGTH * sizeof(uint64_t *));
            for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
                auth_path_a[i] = (uint64_t *)malloc(HASH_SIZE);  // Allocate space for each hash
            }
            size_t proof_len_a = 0, proof_len_b = 0, proof_len_c = 0;
            size_t *proof_len_ptr_a = (size_t *)pull_object(global_proof_stream);
            pull_count++;
            proof_len_a = *proof_len_ptr_a;
            for (size_t i = 0; i < proof_len_a; i++) {
                auth_path_a[i] = (uint64_t *)pull_object(global_proof_stream);
                pull_count++;
            }
            for (size_t i = 0; i < proof_len_a; i++) { //just a debug print
                printf("Layer %zu: ", i);
                if (i < proof_len_a - 1) {  // Codeword + Hash
                    printf("Codeword || Hash: ");
                    print_element(auth_path_a[i], 2 * field_words);  // Concatenated codeword and hash
                } else {  // Only the final hash at the root level
                    printf("Hash: ");
                    print_element(auth_path_a[i], HASH_WORDS);
                }
            }
            printf("pull in verify for auth_path_a %d\n", pull_count);
            if (!merkle_verify(root_verify, aa[s], global_query_indices.a_indices[s], auth_path_a, proof_len_a, field_words)) {
                printf("merkle authentication path verification fails for aa\n");
                free(aa);
                free(bb);
                free(cc);
                free(root_verify);
                for (int j = 0; j < last_domain.size; j++) {
                    free(last_domain_elements[j]);
                }
                free(last_domain_elements);
                return 0;
            }

            //verify auth_path_b
            uint64_t **auth_path_b = (uint64_t **)malloc(MAX_PROOF_PATH_LENGTH * sizeof(uint64_t *));
            for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
                auth_path_b[i] = (uint64_t *)malloc(HASH_SIZE);  // Allocate space for each hash
            }
            size_t *proof_len_ptr_b = (size_t *)pull_object(global_proof_stream);
            proof_len_b = *proof_len_ptr_b;
  
            for (size_t i = 0; i < proof_len_b; i++) {
                auth_path_b[i] = (uint64_t *)pull_object(global_proof_stream); // Pull each hash
            }
            pull_count++;
            printf("pull in verify for auth_path_b %d\n", pull_count);
            if (!merkle_verify(root_verify, bb[s],global_query_indices.b_indices[s], auth_path_b, proof_len_b, field_words)) {
                printf("merkle authentication path verification fails for bb\n");
                free(aa);
                free(bb);
                free(cc);
                //free(roots);
                //free(alphas);
                for (int j = 0; j < last_domain.size; j++) {
                    free(last_domain_elements[j]);
                }
                free(last_domain_elements);
                return 0;
            }

            //verify auth_path_c
            uint64_t **auth_path_c = (uint64_t **)malloc(MAX_PROOF_PATH_LENGTH * sizeof(uint64_t *));
            for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
                auth_path_c[i] = (uint64_t *)malloc(HASH_SIZE);  // Allocate space for each hash
            }
            size_t *proof_len_ptr_c = (size_t *)pull_object(global_proof_stream);
            proof_len_c = *proof_len_ptr_c;

            for (size_t i = 0; i < proof_len_c; i++) {
                auth_path_c[i] = (uint64_t *)pull_object(global_proof_stream); // Pull each hash
            }
            pull_count++;
            printf("pull in verify for auth_path_b %d\n", pull_count);
            if (!merkle_verify(root_verify, cc[s], global_query_indices.c_indices[s], auth_path_c, proof_len_c, field_words)) {
                printf("merkle authentication path verification fails for cc\n");
                free(aa);
                free(bb);
                free(cc);
                //free(roots);
                //free(alphas);
                for (int j = 0; j < last_domain.size; j++) {
                    free(last_domain_elements[j]);
                }
                free(last_domain_elements);
                return 0;
            }
        }
        free(global_query_indices.c_indices);
        free(global_query_indices.a_indices);
        free(global_query_indices.b_indices);
        free(aa);
        free(bb);
        free(cc);

    }

    printf("Verification successful!\n");
    return 1;

}

void test_fri(){
    int degree = 4095;
    int expansion_factor = 32;
    int num_colinearity_tests = 13;

    int initial_domain_length = (degree + 1) * expansion_factor;
    int log_codeword_length = (int)log2(initial_domain_length);
    int basis_len = log_codeword_length;
    uint64_t *eval_basis = populate_eval_basis(basis_len);
    uint64_t offset = get_offset_for_basis(basis_len);
    assert(eval_basis != NULL);
    assert(offset != NULL);

    Fri *fri = init_fri(initial_domain_length, expansion_factor, num_colinearity_tests);
    assert(fri != NULL);

    printf("Number of rounds: %d\n", fri_num_rounds(fri));

    //initialize polynomial coefficients
    int total_elements;
    uint64_t *all_bases = populate_all_bases(&total_elements);

    uint64_t poly_coeffs[initial_domain_length >> fri_num_rounds(fri)];
    for(int i=0; i< initial_domain_length >> fri_num_rounds(fri); i++){
        memcpy(&poly_coeffs[i], &all_bases[i], sizeof(uint64_t));
    }
    uint64_t **domain_elements = (uint64_t **)malloc(initial_domain_length * sizeof(uint64_t *));

    for (int i = 0; i < initial_domain_length; i++) {
        domain_elements[i] = (uint64_t *)malloc(field_words * sizeof(uint64_t));
        assert(domain_elements[i] != NULL);
        uint64_t temp[field_words];
        i_th_ele_in_span(temp, eval_basis, MAX_EVAL_BASIS_LEN, i);
        //print_field("temp", temp, FIELD_WORDS);
        field_add(domain_elements[i], &offset, temp, field_words);
        //print_field("domain_elements[i]", domain_elements[i], FIELD_WORDS);
    }
    // Evaluate polynomial over the domain
    uint64_t **codeword = (uint64_t **)malloc(initial_domain_length * sizeof(uint64_t *));

    for (int i = 0; i < initial_domain_length; i++) {
        codeword[i] = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
        // print_field("codeword value in test_fri", codeword[i], field_words);
        // int is = is_zero(codeword[i], field_words);
    }

    parallel_poly_eval(codeword, poly_coeffs, 32, domain_elements, FIELD_WORDS, initial_domain_length);
    printf("Codeword length: %d\n", initial_domain_length);

    // Test valid codeword
    printf("Testing valid codeword ...\n");
    printf("Calling prove\n");

    prove(fri, codeword);
    printf("Returned from prove for valid codeword\n");

    //size_t num_points = 0;
    //what are these? this is the degree of the polynomial that professor was talking about?
    uint64_t *points[32768];
    for (int i = 0; i < 32768; i++) {
        points[i] = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
    }
    int verdict = verify(fri, points, degree);

    if (verdict == 1)
        printf("accept proof\n");
    else {
        //printf("rejecting proof, but proof should be valid!");
        return;
    }
    
    uint64_t **points_x = (uint64_t **)malloc(13 * sizeof(uint64_t *));
    uint64_t **points_y = (uint64_t **)malloc(13 * sizeof(uint64_t *));
    for (int i = 0; i < 13; i++) {
        points_x[i] = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
        points_y[i] = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
        assert(points_x[i] != NULL);
        assert(points_y[i] != NULL);

        // Assuming points array is structured such that points_x[i] = points[2*i] and points_y[i] = points[2*i+1]
        memcpy(points_x[i], points[2 * i], FIELD_WORDS * sizeof(uint64_t));
        memcpy(points_y[i], points[2 * i + 1], FIELD_WORDS * sizeof(uint64_t));
    }
    int max_basis_len = sizeof(preon.fri_domains) / sizeof(preon.fri_domains[0]);
    for (int i = 0; i < 13; i++) {
        uint64_t temp[FIELD_WORDS];
        i_th_ele_in_span(temp, eval_basis, MAX_EVAL_BASIS_LEN, i);
        field_add(temp, &offset, temp, field_words);

        uint64_t eval_result[FIELD_WORDS];
        poly_eval(eval_result, poly_coeffs, max_basis_len, temp, field_words);

        if (!field_equal(eval_result, points_y[i], FIELD_WORDS)) {
            printf("Polynomial evaluates to wrong value\n");
            assert(false);
        }
    }

    printf("Success! \\o/\n");

    // Disturb then test for failure
    printf("Testing invalid codeword ...\n");
    for (int i = 0; i < degree / 3; i++) {
        memset(codeword[i], 0, FIELD_WORDS * sizeof(uint64_t));
    }

    printf("Calling prove for invalid codeword\n");
    prove(fri, codeword);
    printf("Returned from prove for invalid codeword\n");

    if (verify(fri, points, degree) != 0) {
        fprintf(stderr, "Proof should fail, but is accepted ...\n");
        exit(EXIT_FAILURE);
    }
    printf("Success! \\o/\n");

}

int main() {
    int basis_len = 20;
    uint64_t *eval_basis = populate_eval_basis(basis_len);

    if (eval_basis != NULL) {
        printf("Eval basis for basis_len %d:\n", basis_len);
        for (int i = 0; i < basis_len; i++) {
            printf("%016lx ", eval_basis[i]);
        }
        printf("\n");

        free(eval_basis); 
    }

    uint64_t offset = get_offset_for_basis(basis_len);

    if (offset != 0) {
        printf("Offset for basis_len %d is: %016lx\n", basis_len, offset);
    } else {
        printf("Offset not found or is zero.\n");
    }

    int total_elements;
    uint64_t *all_bases = populate_all_bases(&total_elements);

    if (all_bases != NULL) {
        printf("Total non-zero basis elements: %d\n", total_elements);
        for (int i = 0; i < total_elements; i++) {
            printf("%016lx ", all_bases[i]);
        }
        printf("\n");

        //free all_bases
        free(all_bases);
    } else {
        printf("Error: Could not populate basis elements.\n");
    }

    clock_t t; 
    t = clock(); 
    test_fri(); 
    t = clock() - t; 
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
 
    printf("fri took %f seconds to execute \n", time_taken); 

    return 0;
}

void int_to_bytes(int n, unsigned char *bytes, size_t size) {
    for (size_t i = 0; i < size; i++) {
        bytes[size - 1 - i] = (unsigned char)(n >> (i * 8));
    }
}
