#ifndef COMMIT_KERNEL_CUH
#define COMMIT_KERNEL_CUH

#include <vector>
#include <libff/algebra/fields/binary/gf256.hpp>

using FieldT = libff::gf256;

// Launch FRI commit kernel on GPU (Cantor basis, denominator-free)
// Computes folded codeword: next[j] = current[2j] + alpha * (current[2j+1] - current[2j])
// Parameters:
//   codeword:         current codeword (input, size N)
//   codeword_nxt:     next codeword (output, size N/2)
//   alpha:            challenge field element
//   N:                current codeword size (must be even)
//   layer_num:        merkle tree layer number (0-20)
//   merkle_tree_layers: merkle tree storage (all layers)
void commit_launch(
    uint64_t** codeword,
    uint64_t** codeword_nxt,
    uint64_t* alpha,
    int N,
    uint64_t *root,
    uint64_t** tree_layer,
    uint64_t** tree_layer_nxt,
    uint64_t*** tree,
    int last_round,
    bool is_last_round
);

#endif // COMMIT_KERNEL_CUH
