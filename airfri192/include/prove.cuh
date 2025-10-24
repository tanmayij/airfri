#include "iostream"
#include "vector"
#include "fri_utils.cuh"
using namespace std;

//function definitions 
#include <vector>
#include <libff/algebra/fields/binary/gf256.hpp>
using FieldT = libff::gf256;
uint64_t ***commit_host(Fri *fri, uint64_t **codeword, int codeword_len, uint64_t **tree_layer);
size_t* prove(Fri* fri, uint64_t **codeword, uint64_t **tree_layer);
