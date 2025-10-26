
#include "../include/prove.cuh"
#include "fs_transform.cuh"
#include <openssl/evp.h>
using namespace std;
ProofStream* global_proof_stream = NULL;
#define FIELD_WORDS 4 //gf256 uses 4 uint64_t words
#include "commit_kernel.cuh"
#include "fs_transform.cuh"
#include "sample_indices.cuh"
#include "fri_utils.cuh"
#include "../include/hash_host.cuh"
#include "../include/query.cuh"
#include <iostream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include "../include/merkle.cuh"
#include <libff/algebra/fields/binary/gf256.hpp>
#include <libiop/algebra/field_subset/field_subset.hpp>

#define HASH_WORDS 4

using namespace std;
typedef libff::gf256 FieldT;


extern ProofStream* global_proof_stream;
void int_to_bytes(int n, unsigned char *bytes, size_t size) {
    for (size_t i = 0; i < size; i++) {
        bytes[size - 1 - i] = (unsigned char)(n >> (i * 8));
    }
}
// void hash_sha3_256(const uint64_t *data, size_t len, uint64_t *out) {
//     EVP_MD_CTX *mdctx = EVP_MD_CTX_new();
//     const EVP_MD *md = EVP_sha3_256();

//     EVP_DigestInit_ex(mdctx, md, NULL);
//     EVP_DigestUpdate(mdctx, data, len);
//     EVP_DigestFinal_ex(mdctx, (unsigned char *)out, NULL);

//     EVP_MD_CTX_free(mdctx);
// }
uint64_t ***commit_host(Fri *fri, uint64_t **codeword, int codeword_len, uint64_t **tree_layer) {
    cout << "=== Commit Phase (Host) ===" << endl;
    // uint64_t ***codewords = (uint64_t ***)malloc(fri_num_rounds(fri) + 1 * sizeof(uint64_t **));
    uint64_t ***codewords = new uint64_t**[15 + 1]; //remove hardcoded val later 
    int last_round = fri_num_rounds(fri);
    //get random bytes from Fiat-Shamir proof stash
    size_t N = codeword_len;
    // size_t num_bytes_alpha = 24; //unused variable, removed to fix warning
    // Sample alpha bytes deterministically from the proof stream
    // unsigned char* alpha_bytes = prover_fiat_shamir(global_proof_stream, num_bytes_alpha);
    //use libiop to sample field element from bytes
    //libiop::field_subset<FieldT> can sample random elements
    FieldT alpha = FieldT::random_element();
    // Convert alpha to 4x uint64_t words using to_words()
    std::vector<uint64_t> alpha_vec = alpha.to_words();
    uint64_t *alpha_words = new uint64_t[FIELD_WORDS];
    for (int i = 0; i < FIELD_WORDS; ++i) alpha_words[i] = alpha_vec[i];
    // push_object(global_proof_stream, (void*)alpha_words);
    //allocate space for root hash (will be computed on GPU)
    uint64_t *root = new uint64_t[HASH_WORDS];
    uint64_t *alpha_hash_output = new uint64_t[HASH_WORDS];
    memset(root, 0, HASH_WORDS * 8);
    memset(alpha_hash_output, 0, 32);
    cout << "Alpha sampled and root hash space allocated." << endl;
    cout << "Codeword computations will happen on GPU kernel." << endl;
    for(int r = 0; r < fri_num_rounds(fri); r++) {
        cout << "Round " << r << " processing..." << endl;
        hash_sha3_256(alpha_words, 32, alpha_hash_output);
        // Optionally update alpha_words for next round if needed (not shown)
        cout << "Alpha hashed for round " << r << "." << endl;
        for(int i = 0; i < HASH_WORDS; i++) {
            cout << "Alpha hash output word " << i << ": " << ((uint64_t*)alpha_hash_output)[i] << endl;
        }
        bool is_last_round = (r == fri_num_rounds(fri) - 1);
        cout << "Is last round: " << is_last_round << endl;
        codewords[r] = codeword;
        tree[r] = tree_layer;
        uint64_t **codeword_nxt = new uint64_t*[N/2];
        uint64_t **tree_layer_nxt = new uint64_t*[N/2];
        for(int i = 0; i < N/2; i++) {
            codeword_nxt[i] = new uint64_t[FIELD_WORDS];
            tree_layer_nxt[i] = new uint64_t[FIELD_WORDS];
        }
        // Launch GPU kernel to compute next codeword and tree layer
        // void commit_launch(
        //     uint64_t** codeword,
        //     uint64_t** codeword_nxt,
        //     uint64_t* alpha,
        //     int N,
        //     int layer_num,
        //     uint64_t** merkle_tree_layers
        // );
        commit_launch(codeword, codeword_nxt, alpha_words, N, root, tree_layer, tree_layer_nxt, tree, last_round, is_last_round);
        if(N > 32) {
            N = N / 2;
            codeword = codeword_nxt;
            tree_layer = tree_layer_nxt;
            codeword_len = codeword_len /2;
        }
    }
    
    // After all rounds complete, push to proof stream (ONCE)
    cout << "Pushing to proof stream:" << endl;
    
    // 1. Push merkle root
    push_object(global_proof_stream, (void*)root);
    push_count++;
    cout << "  Pushed root hash" << endl;
    
    // 2. Push alpha
    push_object(global_proof_stream, (void*)alpha_words);
    push_count++;
    cout << "  Pushed alpha" << endl;
    
    // 3. Push last codeword (final codeword after all folding)
    // Final codeword length should be 32 (2^5)
    int final_codeword_length = 32;
    cout << "  Pushing last codeword (" << final_codeword_length << " elements)" << endl;
    for (int i = 0; i < final_codeword_length; i++) {
        push_object(global_proof_stream, (void*)codeword[i]);
        push_count++;
    }
    // cout << "  Pushed " << final_codeword_length << " last codeword elements" << endl;
    
    cout << "Total push count in commit: " << push_count << endl;
    
    return codewords;
}

size_t* prove(Fri* fri, uint64_t **codeword, uint64_t **tree_layer) {
    // You must pass the codeword length and tree layer length as needed
    // For now, assume initial domain length for codeword
    size_t codeword_length = fri->initial_domain_length;
    cout << "codeword length in prove(): " << codeword_length << endl;

    // Create proof stream and run commit phase
    global_proof_stream = create_proof_stream();
    uint64_t ***codewords = commit_host(fri, codeword, codeword_length, tree_layer);

    // Sample indices using Fiat-Shamir (deterministic from transcript)
    size_t num_bytes = 32; // SHA3-256 hash size
    uint8_t *seed = prover_fiat_shamir(global_proof_stream, num_bytes);

    // Allocate space for indices
    size_t *top_level_indices = new size_t[fri->num_colinearity_tests];
    size_t *reduced_indices = new size_t[fri->num_colinearity_tests];

    // Sample indices deterministically from seed
    sample_indices(seed, 32,
                   codeword_length / 2,  // size (after first fold)
                   fri->initial_domain_length >> (fri_num_rounds(fri) - 1),  // reduced_size
                   fri->num_colinearity_tests,  // number of indices
                   top_level_indices,  // output: indices
                   reduced_indices);   // output: reduced indices

    cout << "Sampled " << fri->num_colinearity_tests << " indices from Fiat-Shamir:" << endl;
    for (size_t i = 0; i < fri->num_colinearity_tests; i++) {
        cout << "  index[" << i << "]: " << top_level_indices[i] << endl;
    }

    // Copy indices for query phase
    size_t *indices = new size_t[fri->num_colinearity_tests];
    memcpy(indices, top_level_indices, fri->num_colinearity_tests * sizeof(size_t));

    // Run query phase for each FRI round-- optimize this. it is enough to run first round. 
    //for (size_t r = 0; r < fri_num_rounds(fri) - 1; r++) {
        // Query this round with current indices
    indices = query(fri, codewords,
                codewords[0], codeword_length >> 0,  // current codeword
                codewords[1], codeword_length >> (1),  // next codeword
                indices, 0);
    //}
    cout << "Completed query phase for all rounds." << endl;
    free(seed);
    free(reduced_indices);

    return top_level_indices;
}