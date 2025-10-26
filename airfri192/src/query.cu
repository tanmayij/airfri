#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "../include/query.cuh"
#include "../include/fri_utils.cuh"
#include "../include/fs_transform.cuh"
#include "../include/merkle.cuh"
#include "../include/prove.cuh"
using namespace std;
QueryIndices global_query_indices;
void calculate_indices(size_t  *c_indices, size_t  *a_indices, size_t  *b_indices, int num_colinearity_tests) {

    for (int i = 0; i < num_colinearity_tests; i++) {
        a_indices[i] = (2 * c_indices[i]);
        b_indices[i] = (2 * c_indices[i] + 1);
    }
}

size_t* query(Fri *fri, uint64_t ***codewords, uint64_t **current_codeword, size_t current_codeword_len, uint64_t **next_codeword, size_t next_codeword_len, size_t *c_indices, int round) 
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
        // Print as FieldT vector (assumes current_codeword stores uint64_t[4] for GF(2^256))
        {
            std::vector<uint64_t> words(current_codeword[global_query_indices.a_indices[s]], current_codeword[global_query_indices.a_indices[s]] + field_words);
            FieldT elem(words[0], words[1], words[2], words[3]);
            std::cout << "current_codeword[a_indices[s]]: " << elem << std::endl;
        }
        push_count++;
        push_object(global_proof_stream, current_codeword[global_query_indices.b_indices[s]]);
        {
            std::vector<uint64_t> words(current_codeword[global_query_indices.b_indices[s]], current_codeword[global_query_indices.b_indices[s]] + field_words);
            FieldT elem(words[0], words[1], words[2], words[3]);
            std::cout << "current_codeword[b_indices[s]]: " << elem << std::endl;
        }
        push_count++;
        push_object(global_proof_stream, next_codeword[global_query_indices.c_indices[s]]);
        {
            std::vector<uint64_t> words(next_codeword[global_query_indices.c_indices[s]], next_codeword[global_query_indices.c_indices[s]] + field_words);
            FieldT elem(words[0], words[1], words[2], words[3]);
            std::cout << "next_codeword[c_indices[s]]: " << elem << std::endl;
        }
        push_count++;

        //below logic is to compute auth_paths for each of the codewords and their indices.

        //for a_index
    size_t proof_len_a = 0;
    uint64_t **auth_path_a = new uint64_t*[MAX_PROOF_PATH_LENGTH];
    size_t *proof_len_ptr_a = new size_t;
        // for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) { 
        // auth_path_a[i] = (uint64_t *)malloc(CONCAT_WORDS * sizeof(uint64_t));  
        // }
        merkle_open(auth_path_a, global_query_indices.a_indices[s], &proof_len_a, tree, round);
        *proof_len_ptr_a = proof_len_a;
        push_object(global_proof_stream, proof_len_ptr_a);
        push_count++;
        if (auth_path_a[0]) {
            printf("auth_path_a[0]: %016lx ", *auth_path_a[0]);
        } else {
            printf("auth_path_a[0]: NULL ");
        }
        for (size_t i = 0; i < proof_len_a; i++) {
        //auth_path_a[%zu]\n", i);
        push_object(global_proof_stream, auth_path_a[i]); // Push each hash
        push_count++;
        }
        printf("proof len a: %zu \n", proof_len_a);
        printf("push counted +1 in query for a_indices %d\n", push_count);
        // for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
        //     free(auth_path_a[i]);
        // }
    delete[] auth_path_a;
    // delete proof_len_ptr_a;

        //for b_index
    size_t proof_len_b = 0; 
    uint64_t **auth_path_b = new uint64_t*[MAX_PROOF_PATH_LENGTH];
    size_t *proof_len_ptr_b = new size_t;
        // for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
        // auth_path_b[i] = (uint64_t *)malloc(HASH_SIZE);  //allocate space for each hash
        // }
        merkle_open(auth_path_b, global_query_indices.b_indices[s], &proof_len_b, tree, round);
        *proof_len_ptr_b = proof_len_b;
        push_object(global_proof_stream, proof_len_ptr_b);
               push_count++;
        for (size_t i = 0; i < proof_len_b; i++) {
        //printf("Pushing auth_path_b[%zu]\n", i);
        push_object(global_proof_stream, auth_path_b[i]); //push each hash
        push_count++;
        }
        printf("proof len b: %zu \n", proof_len_b);
        printf("push counted +1 in query for b_indices %d\n", push_count);
        // for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
        //     free(auth_path_b[i]);
        // }
    delete[] auth_path_b;
        

        //for c_index
    size_t proof_len_c = 0;
    uint64_t **auth_path_c = new uint64_t*[MAX_PROOF_PATH_LENGTH];
    size_t *proof_len_ptr_c = new size_t;
        // for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
        // auth_path_c[i] = (uint64_t *)malloc(HASH_SIZE);  //allocate space for each hash
        // }
        merkle_open(auth_path_c, global_query_indices.c_indices[s], &proof_len_c, tree, round+1);
        *proof_len_ptr_c = proof_len_c;
        push_object(global_proof_stream, proof_len_ptr_c);
        push_count++;
        for (size_t i = 0; i < proof_len_c; i++) {
        //printf("Pushing auth_path_c[%zu]\n", i);
        push_object(global_proof_stream, auth_path_c[i]); //push each hash
        push_count++;
        }
        printf("proof len c: %zu \n", proof_len_c);
        // for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH - 1; i++) {
        //     free(auth_path_c[i]);
        // }
    delete[] auth_path_c;
        // free(proof_len_ptr_c);
    }
    printf("Indices in query:\n");
    for (int i = 0; i < fri->num_colinearity_tests; i++) {
    printf("a: %zu, b: %zu, c: %zu\n", global_query_indices.a_indices[i], global_query_indices.b_indices[i], global_query_indices.c_indices[i]);
    }
    return c_indices;
}