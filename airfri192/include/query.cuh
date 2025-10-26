//add headers
#include "prove.cuh"
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "../include/fri_utils.cuh"

void calculate_indices(size_t  *c_indices, size_t  *a_indices, size_t  *b_indices, int num_colinearity_tests);
size_t* query(Fri *fri, uint64_t ***codewords, uint64_t **current_codeword, size_t current_codeword_len, uint64_t **next_codeword, size_t next_codeword_len, size_t *c_indices, int round);
