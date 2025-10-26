#include "fri_utils.cuh"
#include <cmath>
#include <cstdint>
size_t num_leaves = 1048576; // 2^20
int push_count = 0;
int pull_count = 0;
int field_words = 4;
uint64_t ***tree = new uint64_t**[MAX_PROOF_PATH_LENGTH + 1];

int fri_num_rounds(Fri* fri) {
    int codeword_length = fri->initial_domain_length;
    int num_rounds = 0;
    int max_rounds = 15;
    while (num_rounds < max_rounds && codeword_length > fri->expansion_factor && fri->num_colinearity_tests < codeword_length) {
        codeword_length /= 2;
        num_rounds += 1;
    }
    return num_rounds;
}
