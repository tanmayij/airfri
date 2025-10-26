#include <cmath>
#include <cstdint>
#ifndef FRI_UTILS_CUH
#define FRI_UTILS_CUH
extern size_t num_leaves;
#define NUM_LEAVES num_leaves
#define MAX_PROOF_PATH_LENGTH ((int)log2(NUM_LEAVES)) //20
extern int push_count;
extern int pull_count;
extern int field_words;
extern uint64_t ***tree;
struct Fri {
    int initial_domain_length;
    int expansion_factor;
    int num_colinearity_tests;
    int num_rounds;
};
typedef struct {
    size_t *a_indices;
    size_t *b_indices;
    size_t *c_indices;
    int num_colinearity_tests;
} QueryIndices;

inline Fri* init_fri(int initial_domain_length, int expansion_factor, int num_colinearity_tests) {
    Fri* fri = new Fri;
    fri->initial_domain_length = initial_domain_length;
    fri->expansion_factor = expansion_factor;
    fri->num_colinearity_tests = num_colinearity_tests;
    fri->num_rounds = static_cast<int>(log2(initial_domain_length)) - 5;
    return fri;
}
int fri_num_rounds(Fri* fri);

#endif // FRI_UTILS_CUH
