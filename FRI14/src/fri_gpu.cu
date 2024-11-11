#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <cuda_runtime.h>
#include <math.h>
#include "../include/fiat-shamir.cuh"
#include "../include/fri-gpu.cuh"
#include "../include/merkle.cuh"
#include "../include/params.cuh"
#include "../include/field.cuh"
#include "../include/poly.cuh"
#include "../include/domain.cuh"
#include "../include/commit-launch-merkle.cuh"
#include "../include/poly-eval-launch.cuh"

ProofStream* global_proof_stream = NULL;
QueryIndices global_query_indices;

// Extract FRI parameters excluding zeroes
void extract_fri_parameters(const Domain *domains, int num_domains, uint64_t parameters[MAX_FRI_PARAMETERS][2][MAX_EVAL_BASIS_LEN])
{
    for (int i = 0; i < num_domains; i++)
    {
        const Domain *domain = &domains[i];
        int basis_index = 0;
        for (int j = 0; j < domain->basis_len; j++)
        {
            // Assuming the non-zero elements are stored in the last index for each 256-bit element
            if (domain->basis[j * 4] != 0 || domain->basis[j * 4 + 1] != 0 || domain->basis[j * 4 + 2] != 0 || domain->basis[j * 4 + 3] != 0)
            {
                parameters[i][0][basis_index] = domain->basis[j * 4 + 3];
                basis_index++;
            }
        }
        parameters[i][1][0] = domain->shift[3];  // Assuming the 3rd element is non-zero for the shift
    }
}

// Extract FRI parameters excluding zeroes
void init_fri_parameters(uint64_t parameters[MAX_FRI_PARAMETERS][2][MAX_EVAL_BASIS_LEN])
{
    const Domain *domains = preon.fri_domains;
    int num_domains = sizeof(preon.fri_domains) / sizeof(preon.fri_domains[0]);
    extract_fri_parameters(domains, num_domains, parameters);
}

__host__ __device__ void i_th_ele_in_span(uint64_t *result, uint64_t *basis, int len_basis, int i)
{
    *result = 0;  // Initialize result

    // #ifdef __CUDA_ARCH__
    //     printf("GPU: Input index i: %d, len_basis: %d\n", i, len_basis);
    //     // Print all basis elements
    //     for (int j = 0; j < len_basis; ++j) {
    //         printf("GPU: basis[%d] = %016llx\n", j, (unsigned long long)basis[j]);
    //     }
    // #else
    //     printf("CPU: Input index i: %d, len_basis: %d\n", i, len_basis);
    //     // Print all basis elements
    //     for (int j = 0; j < len_basis; ++j) {
    //         printf("CPU: basis[%d] = %016llx\n", j, (unsigned long long)basis[j]);
    //     }
    // #endif

    // Iterate over each bit of the integer i
    for (int bit = 0; bit < len_basis; ++bit)
    {
        // Check if the bit at position 'bit' is set in the binary representation of i
        if ((i >> bit) & 1)
        {
            // #ifdef __CUDA_ARCH__
            //     printf("GPU: Bit %d is set. Adding basis[%d] = %016llx to result\n", bit, bit, (unsigned long long)basis[bit]);
            // #else
            //     printf("CPU: Bit %d is set. Adding basis[%d] = %016llx to result\n", bit, bit, (unsigned long long)basis[bit]);
            // #endif
            // Add the corresponding basis element to the result (XOR operation for uint64_t element)
            *result ^= basis[bit];  // XOR for single uint64_t element
        }
        // else
        // {
        //     #ifdef __CUDA_ARCH__
        //         printf("GPU: Bit %d is not set.\n", bit);
        //     #else
        //         printf("CPU: Bit %d is not set.\n", bit);
        //     #endif
        // }
    }

    // #ifdef __CUDA_ARCH__
    //     printf("GPU: Final result = %016llx\n", (unsigned long long)*result);
    // #else
    //     printf("CPU: Final result = %016llx\n", (unsigned long long)*result);
    // #endif
}
//also need to populate eval_basis and offset_dict here
void populate_eval_basis_and_offset_dict() {
    // field_words = preon.field_words;
    
    // Initialize FRI parameters
    uint64_t FRI_PARAMETERS[MAX_FRI_PARAMETERS][2][MAX_EVAL_BASIS_LEN] = {0};
    init_fri_parameters(FRI_PARAMETERS);

    // Populate GF_2_192_BASIS with actual basis elements from preon.fri_domains
    uint64_t basis_index = 0;
    int fri_domain_count = sizeof(preon.fri_domains) / sizeof(preon.fri_domains[0]);
    printf("fri_domain_count:%d /*replaced from uint64_t */", fri_domain_count);
    for (int i = 0; i < fri_domain_count; i++) {
    for (int j = 0; j < preon.fri_domains[i].basis_len; j++) {
        if (preon.fri_domains[i].basis[j] != 0ull) {
            GF_2_192_BASIS[basis_index++] = preon.fri_domains[i].basis[j];
        }
    }
    }

    // Compute EVAL_BASIS_DICT and OFFSET_DICT
    // Populate EVAL_BASIS_DICT directly from FRI_PARAMETERS
    for (uint64_t exponent = 0; exponent < MAX_FRI_PARAMETERS; exponent++) {
        int eval_len = 0;
        for (int j = 0; j < MAX_EVAL_BASIS_LEN; j++) {
            if (FRI_PARAMETERS[exponent][0][j] != 0) {
                eval_len++;
                EVAL_BASIS_DICT[exponent][j] = (uint64_t *)malloc(field_words * sizeof(uint64_t));  // Allocate memory
                // Directly copy the basis value to EVAL_BASIS_DICT
                EVAL_BASIS_DICT[exponent][j][0] = FRI_PARAMETERS[exponent][0][j]; // Directly assign the value
            } else {
                break;
            }
        }
    

    // printf("EVAL_BASIS_DICT[%d] contents:\n", exponent);
    // for (int /*replaced from uint32_t */ j = 0; j < eval_len; j++) {
    // printf("Element %d /*replaced from uint32_t */: ", j);
    // for (int /*replaced from uint32_t */ k = 0; k < field_words / 4; k++) {
    //     printf("%d ", EVAL_BASIS_DICT[exponent][j][k]);
    // }
    // printf("\n");
    // }
    //Populate OFFSET_DICT
    OFFSET_DICT[exponent] = (uint64_t *)malloc(field_words * sizeof(uint64_t));
    OFFSET_DICT[exponent][0] = FRI_PARAMETERS[exponent][1][0];  // Directly copy the offset value

    // Process the offset element using i_th_ele_in_span or another appropriate method
    //i_th_ele_in_span(OFFSET_DICT[exponent], (const uint64_t *)GF_2_192_BASIS, 64, FRI_PARAMETERS[exponent][1][0], field_words);
    //Populate OFFSET_DICT using the shift parameter
    //OFFSET_DICT[exponent] = (uint64_t *)malloc(field_words * sizeof(uint64_t));
    //i_th_ele_in_span(OFFSET_DICT[exponent], GF_2_192_BASIS, basis_index, preon.fri_domains[exponent].shift[3], field_words); 
    }
    for (uint64_t exponent = 0; exponent < MAX_FRI_PARAMETERS; exponent++) {
    // Print the populated OFFSET_DICT for verification
    printf("Offset for domain %lu: ", exponent);
    for (size_t j = 0; j < 1; j++) {
        printf("%016lx ", OFFSET_DICT[exponent][j]);
    }
    printf("\n");
    }
    // // Print FRI_PARAMETERS
    // for (int i = 0; i < MAX_FRI_PARAMETERS; i++) {
    //     printf("FRI_PARAMETERS[%d]:\n", i);
    //     printf("  Basis: ");
    //     for (int j = 0; j < MAX_EVAL_BASIS_LEN; j++) {
    //         printf("%016lu ", (uint64_t)FRI_PARAMETERS[i][0][j]);
    //     }
    //     printf("\n");
    //     printf("  Offset: %016lu\n", (uint64_t)FRI_PARAMETERS[i][1][0]);
    // }
}
void print_eval_basis_dict() {
    printf("EVAL_BASIS_DICT contents:\n");
    for (int i = 0; i < MAX_FRI_PARAMETERS; i++) {
        printf("EVAL_BASIS_DICT[%d]:\n", i);
        for (int j = 0; j < MAX_EVAL_BASIS_LEN; j++) {
            if (EVAL_BASIS_DICT[i][j] == NULL) {
                printf("  NULL\n");
                continue;
            }
            printf("  Element %d: ", j);
            //for (int k = 0; k < 1; k++) {  // Adjust this loop if field_words changes
            printf("%016lx ", EVAL_BASIS_DICT[i][j][0]);
            //}
            printf("\n");
        }
    }
    printf("\n");
}
void cleanup_eval_basis_and_offset_dict() {
    for (uint64_t exponent = 0; exponent < MAX_FRI_PARAMETERS; exponent++) {
        for (size_t j = 0; j < MAX_EVAL_BASIS_LEN; j++) {
            free(EVAL_BASIS_DICT[exponent][j]);
        }
        free(OFFSET_DICT[exponent]);
    }
}
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
        uint64_t var2_cube[field_words];
        field_mul(var2_cube, var2_squared, var2, field_words);
    
        for (int i = 0; i <= 63; i++) {
            field_pow(temp, var1, i, field_words);
            field_mul(temp, temp, var2_cube, field_words);
            memcpy(&elements[idx * field_words], temp, field_words * sizeof(uint64_t));
            idx++;
        }
    }
}
Fri* init_fri(int initial_domain_length, int expansion_factor, int num_colinearity_tests) {
    Fri* fri = (Fri*)malloc(sizeof(Fri));
    fri->initial_domain_length = initial_domain_length;
    fri->expansion_factor = expansion_factor;
    fri->num_colinearity_tests = num_colinearity_tests;
    return fri;
}
// int fri_num_rounds(Fri* fri) {
//     int codeword_length = fri->initial_domain_length;
//     printf("Initial codeword length: %d\n", codeword_length);

//     int num_rounds = 0;
//     int max_rounds = 12;  // Set the maximum number of rounds to 12

//     // Loop until we reach 12 rounds or the codeword length becomes too small
//     while (num_rounds < max_rounds && codeword_length > fri->expansion_factor && fri->num_colinearity_tests < codeword_length) {
//         codeword_length /= 2;
//         num_rounds += 1;
//     }

//     printf("Final codeword length after %d rounds: %d\n", num_rounds, codeword_length);
//     return num_rounds;
// }
int fri_num_rounds(Fri* fri) {
    int codeword_length = fri->initial_domain_length;
    printf("codeword length: %d\n", codeword_length);
    int num_rounds = 0;
    while (codeword_length > fri->expansion_factor && fri->num_colinearity_tests < codeword_length) {
        codeword_length /= 2;
        num_rounds += 1;
    }
    printf("num rounds: %d\n", num_rounds);
    return num_rounds;
}
__host__ __device__ void print_field(const char *label, const uint64_t *field, size_t field_words) {
    printf("%s: ", label);
    for (size_t i = 0; i < field_words; i++) {
        printf("%016lx ", field[i]);
    }
    printf("\n");
}
void byte_array_to_192bit_integer(const unsigned char *byte_array, size_t byte_array_len, uint64_t *output) {
    //is this host or device?

    memset(output, 0, FIELD_WORDS * sizeof(uint64_t));
    uint64_t acc = 0;
    
    // Accumulate the bytes into a single 192-bit integer
    for (size_t i = 0; i < byte_array_len; i++) {
        acc = (acc << 8) | byte_array[i];
        
        // Check if acc exceeds 192 bits and shift appropriately
        if ((i + 1) % (256 / 8) == 0 || i == byte_array_len - 1) {
            output[i / (256 / 8)] = acc & 0xffffffffffffffff;
            acc = 0;
        }
    }
}

void field_copy(uint64_t *dest, const uint64_t *src, size_t field_words) {
    memcpy(dest, src, field_words * sizeof(uint64_t));
}

void field_sample(uint8_t *byte_array, size_t byte_array_len, uint64_t *eval_basis, size_t basis_len, uint64_t *result) {
    uint64_t acc[4] = {0};  // 256 bits (4 * 64 bits)

    printf("Byte array: ");
    for (size_t i = 0; i < byte_array_len; ++i) {
        printf("%02x", byte_array[i]);
    }
    printf("\n");

    // Convert byte_array to a large integer accumulator (acc)
    for (uint64_t i = 0; i < byte_array_len; ++i) {
        // Shift the accumulator left by 8 bits and add the current byte
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

    // Mask to ensure acc is at most 256 bits (though acc is already 256 bits, this is just to be sure)
    acc[0] &= 0xFFFFFFFFFFFFFFFF;
    acc[1] &= 0xFFFFFFFFFFFFFFFF;
    acc[2] &= 0xFFFFFFFFFFFFFFFF;
    acc[3] &= 0xFFFFFFFFFFFFFFFF;

    printf("Accumulator after masking: ");
    for (uint64_t i = 0; i < 4; i++) {
        printf("%016lx ", acc[i]);
    }
    printf("\n");

    // Combine the parts of the accumulator to create a single index
    uint64_t combined_index = acc[0] ^ acc[1] ^ acc[2] ^ acc[3];
    combined_index = combined_index % basis_len;  // Ensure the index is within the range of the basis length
    printf("Combined index: %016lx\n", combined_index);

    // Use the combined index to find the element in the span
    i_th_ele_in_span(result, eval_basis, basis_len, combined_index);

    // Print the result for debugging
    printf("Sampled field element: ");
    for (uint64_t i = 0; i < FIELD_WORDS; i++) {
        printf("%016lx ", result[i]);
    }
    printf("\n");
}


void int_to_bytes(int n, unsigned char *bytes, size_t size) {
    for (size_t i = 0; i < size; i++) {
        bytes[size - 1 - i] = (unsigned char)(n >> (i * 8));
    }
}


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

        counter += 1;  // Increment counter

        // Check for uniqueness
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

//i am not sure if the below function will work in GPU
__host__ __device__ void print_codeword(uint64_t **codeword, size_t length, const size_t field_words) {
    for (size_t i = 0; i < length; i++) {
        printf("codeword[%zu /*replaced from uint32_t */]: ", i);
        for (size_t j = 0; j < field_words; j++) {
            printf("%016lx ", codeword[i][j]);
        }
        printf("\n");
    }
}

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

/* COMMIT PHASE */

uint64_t ***commit_host(Fri *fri, uint64_t **codeword, int codeword_len) //are we allowed to use int /*replaced from uint32_t */ here??
{
    int exponent = (MAX_EVAL_BASIS_LEN - fri_log_domain_length(fri)); //this is 20 - 17 = 3 because of the way the data structure is
    printf("exponent in commit: %d", exponent);
    int basis_len = fri_log_domain_length(fri);
    uint64_t eval_basis[basis_len];
    for (size_t i = 0; i < basis_len; i++) {
        eval_basis[i] = *EVAL_BASIS_DICT[exponent][i];
    } 
    uint64_t *offset = (uint64_t *)malloc(sizeof(uint64_t));
    offset = OFFSET_DICT[exponent];
    //print_eval_basis("eval_basis in commit", eval_basis, basis_len, 1);
    //also not sure if assert will work
    assert(eval_basis != NULL && "eval_basis are not correct");
    //when do i make the below line into cudamalloc??
    uint64_t ***codewords = (uint64_t ***)malloc((fri_num_rounds(fri) - 1) * sizeof(uint64_t **));
    uint64_t precomputed_inverses[MAX_FRI_PARAMETERS]; //another array that needs to be sent to the GPU
    load_precomputed_inverses("fri_inverses_256.txt", precomputed_inverses);

    // Create the lower basis (primitive polynomial of GF(2^64))
    //here, we are creating a basis of 256 elements to sample alphas from. 
    uint64_t var1[FIELD_WORDS] = {0x2}; // Represents var1
    uint64_t var2[FIELD_WORDS] = {0x3}; // Represents var2


    uint64_t elements[256 * FIELD_WORDS];
    //generate elements- this function is not implemented for GPU code yet!
    generate_elements(elements, var1, var2, field_words);


    //alpha needs to be sent to gpu too - for every round, we hash alpha and send the hashed alpha 
    uint64_t *alpha = (uint64_t *) malloc (FIELD_WORDS * sizeof(uint64_t));
    uint64_t *sampled_alpha = (uint64_t *) malloc (FIELD_WORDS * sizeof(uint64_t));

    size_t N = codeword_len;
    printf("num rounds in commit are %d\n", fri_num_rounds(fri));
    uint64_t *root = (uint64_t *)malloc(4 * sizeof(uint64_t)); //32 bytes for merkle tree hash
    //sample alpha here
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
    //iteration for rounds- i feel like for every round, send the codeword information to kernel, let it process it and return to this loop. 
    for(int r = 0; r < fri_num_rounds(fri); r++) {
        //memory allocation here needs to be for both CPU and GPU operations
        //uint64_t *root = (uint64_t *)malloc(4 * sizeof(uint64_t)); //32 bytes for merkle tree hash

        for(int i = 0; i<basis_len - r; i++){
            eval_basis[i] = *EVAL_BASIS_DICT[exponent + r][i];
            }
        ////print_eval_basis("eval basis in commit round", eval_basis, MAX_EVAL_BASIS_LEN - (exponent + r), field_words);
        printf("round: %d\n", r);
        offset = OFFSET_DICT[exponent + r];
        print_field("offset in commit", offset, field_words);
        assert(eval_basis!=NULL && "eval_basis in commit are not correct");
        hash_sha3_256((uint64_t *) sampled_alpha, 24, (uint64_t *) alpha);
        memcpy(sampled_alpha, alpha, 24);//repeat hash - the second time, alpha will be rehashed

        printf("Alpha that is hashed in commit: ");
        for (size_t i = 0; i < FIELD_WORDS; i++) {
            printf("%016lx ", alpha[i]);
        }
        printf("\n");
        
        codewords[r] = codeword;
        uint64_t *denominator_inv = &precomputed_inverses[exponent + r];
        print_field("denominator inverse", denominator_inv, 1);
        uint64_t **codeword_nxt = (uint64_t **)malloc((N/2)*sizeof(uint64_t *)); 

        commit_launch(codeword, codeword_nxt, alpha, offset, denominator_inv[0], eval_basis, N, root);
        if(N > 32) {
        N = N / 2;
        codeword = codeword_nxt;
        // Debug print of codeword (before calling merkle_commit)
        codeword_len = codeword_len /2;
        }
    }
    p 
}

/* QUERY PHASE */

void calculate_indices(size_t  *c_indices, size_t  *a_indices, size_t  *b_indices, int num_colinearity_tests) {

    for (int i = 0; i < num_colinearity_tests; i++) {
        a_indices[i] = (2 * c_indices[i]);
        b_indices[i] = (2 * c_indices[i] + 1);
    }
}
size_t* query(Fri *fri, uint64_t ***codewords, uint64_t **current_codeword, size_t current_codeword_len, uint64_t **next_codeword, size_t next_codeword_len, size_t *c_indices) {
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

    // Push queried elements and their authentication paths into the proof stream
for (int s = 0; s < num_tests; s++) {
    push_object(global_proof_stream, current_codeword[global_query_indices.a_indices[s]]);
    print_field("current_codeword[a_indices[s]]", current_codeword[global_query_indices.a_indices[s]], field_words);
    push_count++;
    push_object(global_proof_stream, current_codeword[global_query_indices.b_indices[s]]);
    print_field("current_codeword[b_indices[s]]", current_codeword[global_query_indices.b_indices[s]], field_words);
    push_count++;
    push_object(global_proof_stream, next_codeword[global_query_indices.c_indices[s]]);
    print_field("next_codeword[c_indices[s]]", next_codeword[global_query_indices.c_indices[s]], field_words);
    push_count++;

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
    //print_field("next codeword[0]", &next_codeword[0], field_words);
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

    // Initialize space for top-level and reduced indices
    size_t *top_level_indices = (size_t *)malloc(fri->num_colinearity_tests * sizeof(size_t));
    size_t *reduced_indices = (size_t *)malloc(fri->num_colinearity_tests * sizeof(size_t));

    // Sample indices for the first round
    sample_indices(seed, 32, codeword_length / 2, 
                    (fri->initial_domain_length >> (fri_num_rounds(fri) - 2)),
                    fri->num_colinearity_tests, 
                    top_level_indices, reduced_indices);

    // Store and reuse indices across rounds
    size_t *indices = (size_t *)malloc(fri->num_colinearity_tests * sizeof(size_t));
    memcpy(indices, top_level_indices, fri->num_colinearity_tests * sizeof(size_t));

    for (size_t i = 0; i < fri_num_rounds(fri) - 2; i++) {
        printf("what is taking time?\n");
        query(fri, codewords, codewords[i], codeword_length, codewords[i + 1], codeword_length / 2, indices);

        // Store the sampled indices for this round
        memcpy(global_indices_tracker.round_indices[i], indices, fri->num_colinearity_tests * sizeof(size_t));

        // Fold indices for the next round
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
    int exponent = (MAX_EVAL_BASIS_LEN - fri_log_domain_length(fri)); //this is 20 - 17 = 3 because of the way the data structure is    uint64_t **eval_basis = EVAL_BASIS_DICT[exponent];    printf("exponent: %d\n", exponent);
    printf("exponent: %d\n", exponent);
    int basis_len = fri_log_domain_length(fri);
    uint64_t eval_basis[basis_len];
    for (size_t i = 0; i < basis_len; i++) {
        eval_basis[i] = *EVAL_BASIS_DICT[exponent][i];
    } 
    uint64_t *offset = (uint64_t *)malloc(sizeof(uint64_t));
    offset = OFFSET_DICT[exponent];
    assert(eval_basis != NULL && "eval_basis are not correct!");
    uint64_t *root_verify = (uint64_t *)malloc(HASH_SIZE * sizeof(uint64_t)); 
    size_t last_codeword_length = (int)pow(2, fri_log_domain_length(fri)) >> fri_num_rounds(fri); //should be 32 for FRI192
    printf("last codeword length: %zu\n", last_codeword_length);
    uint64_t **last_codeword_arr = (uint64_t **)malloc(last_codeword_length * sizeof(uint64_t));
    for (int i = 0; i < last_codeword_length; i ++){
        *last_codeword_arr = (uint64_t *)malloc(field_words * sizeof(uint64_t *));
    }   
    uint64_t var1[FIELD_WORDS] = {0x2}; // Represents var1
    uint64_t var2[FIELD_WORDS] = {0x3}; // Represents var2
    uint64_t elements[256*FIELD_WORDS];

    generate_elements(elements, var1, var2, field_words);

    uint64_t *alpha = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
    uint64_t *sampled_alpha = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
    uint8_t *alpha_bytes = verifier_fiat_shamir(global_proof_stream, 24);

    root_verify = (uint64_t *)pull_object(global_proof_stream);
    sampled_alpha = (uint64_t *)pull_object(global_proof_stream);


    pull_count++;
    pull_count++;
    printf("printing FRI Merkle root");
    for(int j = 0; j< FIELD_WORDS; j++){
        printf("%016lx ", root_verify[j]);
    }
    printf("\n");

    //extract last codeword
    for(int i=0; i<last_codeword_length; i++){
        last_codeword_arr[i] = (uint64_t *)pull_object(global_proof_stream);
        pull_count++;
    }
    printf("\n");
    print_codeword(last_codeword_arr, last_codeword_length, field_words);
    printf("fri_num_rounds: %d", fri_num_rounds(fri));
    uint64_t *check_last_root = (uint64_t *)malloc(HASH_WORDS * sizeof(uint64_t));
    uint64_t *copied_last_codeword = (uint64_t *)malloc(HASH_SIZE);
    //memcpy(copied_root, roots[num_r - 1], HASH_SIZE);
    memcpy(copied_last_codeword, last_codeword_arr, HASH_SIZE);

    degree = ((32 * 10) / fri->expansion_factor) - 1; //(len(last_codeword) // self.expansion_factor) - 1 - this is hardcoded rn but fix that  
    printf("Degree: %d\n", degree);
    int r = fri_num_rounds(fri); //12
    uint64_t *last_basis = (uint64_t *)malloc((int)log2(last_codeword_length) * sizeof(uint64_t));
    for (size_t i = 0; i < (int)log2(last_codeword_length); i++) {
    last_basis = EVAL_BASIS_DICT[exponent + r][i];
    }
    uint64_t last_offset = *OFFSET_DICT[exponent + r];
    //assert(last_basis != NULL && "last_basis are not correct!");
    
    Domain last_domain;
    initialize_domain(&last_domain, fri->initial_domain_length >> r, last_basis, &last_offset);
    uint64_t **domain_elements = (uint64_t **)malloc(last_domain.size * sizeof(uint64_t *));
    
    printf("%zu \n", last_domain.size);
    for (int i = 0; i < last_domain.size; i++) {
        domain_elements[i] = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
        i_th_ele_in_span(domain_elements[i], last_basis, last_domain.basis_len, i);
        field_add(domain_elements[i], domain_elements[i], &last_offset, field_words);
    }

    uint64_t poly_coeffs_[FIELD_WORDS * last_domain.size];
    interpolate_domain_single(poly_coeffs_, domain_elements, copied_last_codeword, last_domain.size, FIELD_WORDS);
    print_array("poly_coeffs", poly_coeffs_, last_domain.size * FIELD_WORDS * sizeof(uint64_t));
    uint64_t eval_result[FIELD_WORDS];
    //printf("last domain size:%d /*replaced from uint32_t */\n", last_domain.size);
    for (int i = 0; i < last_domain.size; i++) {
        poly_eval(eval_result, poly_coeffs_, last_domain.size, domain_elements[i], FIELD_WORDS);

        if (memcmp(copied_last_codeword, eval_result, FIELD_WORDS * sizeof(uint64_t)) == 0) {
            printf("Re-evaluated codeword does not match original!\n");
            free(root_verify);
            for (int j = 0; j < last_domain.size; j++) {
                free(domain_elements[j]);
            }
            free(domain_elements);
            return 0;
        }
    }

    printf("Degree checked!!\n");
    uint64_t poly_coeffs[64] = {0};  // Initialize all elements to 0
    poly_coeffs[0] = 1;
    poly_coeffs[1] = 1;
    poly_coeffs[2] = 0;
    poly_coeffs[3] = 1;
    poly_coeffs[63] = 1;  // Set the 64th element (index 63) to 1
    //size_t poly_len = sizeof(poly_coeffs_) / sizeof(poly_coeffs[0]);
    // printf("poly_len: %d /*replaced from uint32_t */\n",poly_len);
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

    for (r = 0; r < fri_num_rounds(fri) - 2; r++) {
        printf("round %d in verify\n", r);
        hash_sha3_256((uint64_t *) sampled_alpha, 24, (uint64_t *) alpha); //hash alpha
        memcpy(sampled_alpha, alpha, 24); //repeat hash
        printf("Alpha that is hashed in verify: ");
        for (uint64_t  i = 0; i < field_words; i++) {
            printf("%016lx ", alpha[i]);
        }
        printf("\n");
        int basis_len = fri_log_domain_length(fri) >> (r);  
        printf("basis len is %d in round: %d", basis_len, r);
        uint64_t eval_basis[basis_len];  //re-define eval_basis for this round

        for (size_t i = 0; i < basis_len; i++) {
            eval_basis[i] = *EVAL_BASIS_DICT[exponent + r][i];
        }

        offset = OFFSET_DICT[exponent + r];  
        assert(eval_basis != NULL && "eval_basis are not correct!");

        size_t *c_indices = global_indices_tracker.round_indices[r];
        calculate_indices(c_indices, 
            global_query_indices.a_indices, 
            global_query_indices.b_indices, 
            fri->num_colinearity_tests);
        
        printf("Indices in verify:\n");
        for (int i = 0; i < fri->num_colinearity_tests; i++) {
            printf("a: %zu, b: %zu, c: %zu\n", global_query_indices.a_indices[i], global_query_indices.b_indices[i], c_indices[i]);
        }
        uint64_t **aa = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **bb = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **cc = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **ay = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **by = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **cy = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        //colinearity check
        uint64_t precomputed_inverses[MAX_FRI_PARAMETERS];
        load_precomputed_inverses("fri_inverses_256.txt", precomputed_inverses);

        uint64_t *denominator_inv = &precomputed_inverses[exponent + r]; 
        for (int s =0; s < fri->num_colinearity_tests; s++) {
            ay[s] = (uint64_t *)pull_object(global_proof_stream);
            pull_count++;
            //printf("pull in verify for ay %d\n", pull_count);
            print_field("ay", ay[s], field_words);
            aa[s] = ay[s];
            by[s] = (uint64_t *)pull_object(global_proof_stream);
            pull_count++;
            //printf("pull in verify for by %d\n", pull_count);
            print_field("by", by[s], field_words);
            bb[s] = by[s];
            cy[s] = (uint64_t *)pull_object(global_proof_stream);
            pull_count++;
            //printf("pull in verify for cy %d\n", pull_count);
            print_field("cy", cy[s], field_words);
            cc[s] = cy[s];
         
            if (r == 0) {
                polynomial_values[global_query_indices.a_indices[s]] = aa[s];
                polynomial_values[global_query_indices.b_indices[s]] = bb[s];
            }
            //MAKE THIS WORK!!!
            // // Printing the polynomial values that are not zero
            // printf("Polynomial values:\n");
            // int k = 0; // Define k outside the loop
            // for (int i = 0; i < fri->num_colinearity_tests; i++) {
            //     // Check if the polynomial value is zero, skip if true
            //     if (is_zero(polynomial_values[i], 1)) {
            //         continue;
            //     }
            
            //     // Print the current non-zero polynomial value
            //     for (int j = 0; j < FIELD_WORDS; j++) {  // Assuming FIELD_WORDS is the length of each polynomial value
            //         printf("%016lx ", polynomial_values[i][j]);
            //     }
            //     printf("\n");
            //     for (int j = 0; j < FIELD_WORDS; j++) {
            //         polynomial_values[k][j] = polynomial_values[i][j];
            //     }
            //     k++;
            // }
            // uint64_t *ax = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
            // uint64_t *bx = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
            // //uint64_t bx[FIELD_WORDS];
            // uint64_t *cx = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
            //uint64_t *ax_temp = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
            //uint64_t *bx_temp = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
            uint64_t ax[FIELD_WORDS] = {0};
            uint64_t ax_temp[FIELD_WORDS] = {0};
            uint64_t bx[FIELD_WORDS] = {0};
            uint64_t bx_temp[FIELD_WORDS] = {0};
            uint64_t alpha_offset[field_words] = {0};
            uint64_t cx[FIELD_WORDS] = {0};
            uint64_t diff1[FIELD_WORDS] = {0};
            uint64_t diff2[FIELD_WORDS] = {0};
            uint64_t temp1[FIELD_WORDS] = {0};
            uint64_t temp2[FIELD_WORDS] = {0};
            //uint64_t *temp2 = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
            uint64_t temp3[FIELD_WORDS] = {0};
            uint64_t temp4[field_words] = {0};
            uint64_t temp5[field_words] = {0};
            uint64_t den_inv[field_words] = {0};
            printf("exponent: %d\n", exponent);
            printf("a_index %zu ", global_query_indices.a_indices[s]);
            //print_eval_basis("eval basis in verify line computation", eval_basis, exponent + r - 1, field_words);
            i_th_ele_in_span(ax_temp, eval_basis,  exponent + r-1, global_query_indices.a_indices[s]);
            printf("ax_temp: %016lx\n", ax_temp[0]);  // For one element
            //print_field("ax_temp", ax_temp, field_words);
            // printf("offset: ");
            print_field("offset in verify", offset, field_words);
            field_add(ax, ax_temp, offset, field_words);
            print_field("ax", ax, field_words);
            i_th_ele_in_span(bx_temp, eval_basis,  exponent + r-1, global_query_indices.b_indices[s]);
            printf("bx_temp: %016lx\n", bx_temp[0]);;  // For one element
            //print_field("bx", bx, field_words);
            field_add(bx, offset, bx_temp, field_words);
            // bx[1] = {0};
            // bx[2] = {0};
            // bx[3] = {0};
            print_field("offset in verify", offset, field_words);

            print_field("bx", bx, field_words);
            memcpy(cx, alpha, FIELD_WORDS * sizeof(uint64_t)); //should this be alpha or sampled_alpha
            print_field("cx", cx, field_words);
            field_sub(diff1, ay[s], by[s], field_words);
            print_field("diff1", diff1, field_words);
            //Compute diff2 = ax - bx
            field_sub(diff2, ax, bx, field_words);
            print_field("diff2", diff2, field_words);
            //Compute temp1 = inverse(diff2)
            field_inv(temp1, diff2, FIELD_WORDS);
            memcpy(alpha_offset, ax, sizeof(uint64_t));
            print_field("temp1", denominator_inv, field_words);
            // Compute temp2 = cx - ax
            print_field("alpha_offset", alpha_offset, field_words);
            field_sub(temp2, cx, alpha_offset, field_words);
            print_field("temp2", temp2, field_words);
            field_mul(temp3, temp2, &denominator_inv[0], field_words);
            field_mul(temp4, temp3, diff1, FIELD_WORDS);
            // Compute intermediate = diff1 * denominator_inv * temp2
            //field_mul(temp3, temp2, denominator_inv, field_words);
            print_field("denominator inverse in verify", denominator_inv, 1);
            print_field("temp3", temp3, field_words);
            // field_mul(temp4, temp2, diff1, field_words);
            print_field("temp4", temp4, field_words);
            // Compute temp3 = intermediate + ay
            field_add(temp5, temp4, ay[s], FIELD_WORDS);
            print_field("temp5", temp5, field_words);
            printf("s!!! %d, r? %d \n",s, r);
            if (memcmp(cy[s], temp5, field_words * sizeof(uint64_t)) != 0) {
                printf("Colinearity check failed\n");
                printf("%d\n", memcmp(cy[s], temp5, field_words * sizeof(uint64_t)));
                // free(aa);
                // free(bb);
                // free(cc);
                // //free(root_verify);
                // //free(alpha);
                // for (int j = 0; j < last_domain.size; j++) {
                //     free(domain_elements[j]);
                // }
                // free(domain_elements);
                return 0;
            }
          
   
            uint64_t **auth_path_a = (uint64_t **)malloc(MAX_PROOF_PATH_LENGTH * sizeof(uint64_t *));
            for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
                auth_path_a[i] = (uint64_t *)malloc(HASH_SIZE);  // Allocate space for each hash
            }
            size_t proof_len_a, proof_len_b, proof_len_c;
            size_t *proof_len_ptr_a = (size_t *)pull_object(global_proof_stream);
            pull_count++;
            proof_len_a = *proof_len_ptr_a;

            for (size_t i = 0; i < proof_len_a; i++) {
                auth_path_a[i] = (uint64_t *)pull_object(global_proof_stream); // Pull each hash
                pull_count++;
            }
            for (size_t i = 0; i < proof_len_a; i++) {
                printf("Layer %zu: ", i);
        
                if (i < proof_len_a - 1) {  // Codeword + Hash
                    printf("Codeword || Hash: ");
                    print_element(auth_path_a[i], 2 * field_words);  // Concatenated codeword and hash
                } else {  // Only the final hash at the root level
                    printf("Hash: ");
                    print_element(auth_path_a[i], HASH_WORDS);
                }
            }
            // pull_count++;
            printf("pull in verify for auth_path_a %d\n", pull_count);
            //merkle_verify(root_verify, aa[s], global_query_indices.a_indices[s], auth_path_a, proof_len, field_words)
            if (!merkle_verify(root_verify, aa[s], global_query_indices.a_indices[s], auth_path_a, proof_len_a, field_words)) {
                printf("merkle authentication path verification fails for aa\n");
                free(aa);
                free(bb);
                free(cc);
                free(root_verify);
                for (int j = 0; j < last_domain.size; j++) {
                    free(domain_elements[j]);
                }
                free(domain_elements);
                return 0;
            }
      
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
                    free(domain_elements[j]);
                }
                free(domain_elements);
                return 0;
            }
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
                    free(domain_elements[j]);
                }
                free(domain_elements);
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


void print_eval_basis(const char *label, uint64_t *last_basis, size_t basis_len, size_t field_words) {
    printf("%s:\n", label);
    for (size_t i = 0; i < basis_len; i++) {
        printf("%016lx ", last_basis[i]);
        printf("\n");
    }
}
void print_array(const char *label, const uint64_t *array, size_t size) {
    printf("%s: ", label);
    for (size_t i = 0; i < size / sizeof(uint64_t); i++) {
        printf("%016lx \n", array[i]);
    }
    printf("\n");
}

void test_fri() {
    int degree = 4095;
    int expansion_factor = 32;
    int num_colinearity_tests = 12;

    int initial_domain_length = (degree + 1) * expansion_factor;// 2^17
    int log_codeword_length = (int)log2(initial_domain_length); //17
    printf("log codeword len %d\n", log_codeword_length);
    // Correctly initialize eval_basis and offset
    int basis_len = log_codeword_length;
    uint64_t eval_basis[log_codeword_length];
    int pos = MAX_EVAL_BASIS_LEN - log_codeword_length;
    for (size_t i = 0; i < basis_len; i++) {
        eval_basis[i] = *EVAL_BASIS_DICT[pos][i];
    } 
    uint64_t *offset = (uint64_t *)malloc(sizeof(uint64_t));
    offset = OFFSET_DICT[pos];
    assert(eval_basis != NULL);
    assert(offset != NULL);
    //print_eval_basis("eval_basis", eval_basis, MAX_EVAL_BASIS_LEN, FIELD_WORDS);
    Fri* fri = init_fri(initial_domain_length, expansion_factor, num_colinearity_tests);
    assert(fri != NULL);

    printf("Number of rounds: %d\n", fri_num_rounds(fri));
    printf("Exponent: %d\n", fri_log_domain_length(fri)); //we start here

    // Initialize polynomial coefficients
    uint64_t poly_coeffs[BASIS_SIZE];
    for (int i = 0; i < BASIS_SIZE; i++) {
        memcpy(&poly_coeffs[i], &GF_2_192_BASIS[i], sizeof(uint64_t));
    }

    // Initialize domain
    uint64_t **domain_elements = (uint64_t **)malloc(initial_domain_length * sizeof(uint64_t *));

    for (int i = 0; i < initial_domain_length; i++) {
        domain_elements[i] = (uint64_t *)malloc(field_words * sizeof(uint64_t));
        assert(domain_elements[i] != NULL);
        uint64_t temp[field_words];
        i_th_ele_in_span(temp, eval_basis, MAX_EVAL_BASIS_LEN, i);
        //print_field("temp", temp, FIELD_WORDS);
        field_add(domain_elements[i], offset, temp, field_words);
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

    for (int i = 0; i < 13; i++) {
        uint64_t temp[FIELD_WORDS];
        i_th_ele_in_span(temp, eval_basis, MAX_EVAL_BASIS_LEN, i);
        field_add(temp, offset, temp, field_words);

        uint64_t eval_result[FIELD_WORDS];
        poly_eval(eval_result, poly_coeffs, BASIS_SIZE, temp, field_words);

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
void cleanup(Fri* fri, uint64_t** domain_elements, uint64_t** codeword, uint64_t** points, int num_elements, int num_points, uint64_t** aa, uint64_t** bb, uint64_t** cc, int* c_indices, int* a_indices, int* b_indices, uint64_t** alphas) {
    // Free domain elements
    for (int i = 0; i < num_elements; i++) {
        free(domain_elements[i]);
    }
    free(domain_elements);

    // Free codewords
    for (int i = 0; i < num_elements; i++) {
        free(codeword[i]);
    }
    free(codeword);

    // Free points
    for (int i = 0; i < num_points; i++) {
        free(points[i]);
    }
    free(points);

    // Free FRI struct
    free(fri);

    // Free arrays used in verify
    if (aa) {
        for (int i = 0; i < num_elements; i++) {
            free(aa[i]);
        }
        free(aa);
    }
    if (bb) {
        for (int i = 0; i < num_elements; i++) {
            free(bb[i]);
        }
        free(bb);
    }
    if (cc) {
        for (int i = 0; i < num_elements; i++) {
            free(cc[i]);
        }
        free(cc);
    }
    cleanup_indices_tracker();


    // Free global proof stream
    free_proof_stream(global_proof_stream);

    // Free eval_basis and offset dictionaries
    cleanup_eval_basis_and_offset_dict();
}


int main() {
   populate_eval_basis_and_offset_dict();
    // Fri* fri = NULL;
    clock_t t; 
    t = clock(); 
    test_fri(); 
    t = clock() - t; 
    double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds 
 
    printf("fri took %f seconds to execute \n", time_taken); 
    //cleanup(fri, domain_elements, codeword, points, fri->initial_domain_length, 10, aa, bb, cc, c_indices, a_indices, b_indices, alphas);
    return 0;
}