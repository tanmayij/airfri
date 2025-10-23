#ifndef COMMIT_KERNEL_HPP
#define COMMIT_KERNEL_HPP

#include <vector>
#include <libff/algebra/fields/binary/gf256.hpp>

using FieldT = libff::gf256;

//host function to launch FRI commit kernel on GPU (simplified for cantor basis)
//computes folded codeword: next[j] = current[2j] + alpha * (current[2j+1] - current[2j])
//cantor basis structure eliminates the need for denominator computation
void commit_launch(
    std::vector<FieldT>& codeword,           //current codeword
    std::vector<FieldT>& codeword_nxt,       //next codeword (output)
    FieldT& alpha,                           //challenge field element
    int N,                                   //current codeword size
    int layer_num,                           //merkle tree layer number (0-20)
    std::vector<std::vector<FieldT>>& merkle_tree_layers  //merkle tree storage
);

#endif //COMMIT_KERNEL_HPP
