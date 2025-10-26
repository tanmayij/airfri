#ifndef SAMPLE_INDICES_HPP
#define SAMPLE_INDICES_HPP

#include <cstdint>
#include <cstddef>

/**
 * Sample indices deterministically from a seed using SHA3
 * This ensures verifier can reproduce the same indices
 * 
 * @param seed: Random bytes from Fiat-Shamir
 * @param seed_len: Length of seed (typically 32 bytes)
 * @param size: Size of the domain to sample from (after first fold)
 * @param reduced_size: Size after all folding rounds
 * @param number: Number of indices to sample (num_colinearity_tests)
 * @param indices: Output array for sampled indices
 * @param reduced_indices: Output array for reduced indices
 */
void sample_indices(
    uint8_t* seed, 
    size_t seed_len, 
    int size, 
    int reduced_size, 
    int number, 
    size_t* indices, 
    size_t* reduced_indices
);

#endif // SAMPLE_INDICES_HPP
