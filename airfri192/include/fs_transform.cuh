#pragma once

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <openssl/evp.h>

// Define a structure for the ProofStream
struct ProofStream {
    std::vector<void*> objects;
    size_t read_index;

    ProofStream();
    ~ProofStream();
};

// Initialize a ProofStream
ProofStream* create_proof_stream();

// Push an object into the ProofStream
void push_object(ProofStream *ps, void *obj);

// Pull an object from the ProofStream
void* pull_object(ProofStream *ps);

// Helper function to serialize the ProofStream (basic version)
unsigned char* serialize_proof_stream(ProofStream *ps, size_t num_bytes, size_t *out_size);

// Helper function to compute SHAKE-256
unsigned char* compute_shake256(unsigned char *input, size_t input_len, size_t num_bytes);

// Prover's Fiat-Shamir transformation
unsigned char* prover_fiat_shamir(ProofStream *ps, size_t num_bytes);

// Verifier's Fiat-Shamir transformation
unsigned char* verifier_fiat_shamir(ProofStream *ps, size_t num_bytes);

// Free the proof stream
void free_proof_stream(ProofStream *ps);

// Clear global serialized data
void clear_global_serialized_data();
