#ifndef VERIFY_CUH
#define VERIFY_CUH

#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <cstring>
#include "fri_utils.cuh"
#include "commit_kernel.cuh"
#include "poly_eval.cuh"
#include "merkle.cuh"
#include "fs_transform.cuh"
#include "prove.cuh"

// Forward declarations for types
struct Fri;
struct Domain;
extern ProofStream* global_proof_stream;
extern int pull_count;

// Function declarations
int verify(Fri *fri, uint64_t **polynomial_values, int degree);

#endif // VERIFY_CUH
