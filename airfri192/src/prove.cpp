#include "iostream"
#include "vector"
#include "prove.hpp"
#include "fs_transform.hpp"
using namespace std;
ProofStream* global_proof_stream = NULL;
FIELD_WORDS = 4; //gf256 uses 4 uint64_t words
#include "commit.hpp"
#include "fs_transform.hpp"
#include "sample_indices.hpp"
#include <iostream>
#include <vector>
#include <cstring>
#include <libff/algebra/fields/binary/gf256.hpp>
#include <libiop/algebra/field_subset/field_subset.hpp>

#define HASH_WORDS 4

using namespace std;
typedef libff::gf256 FieldT;

// External global proof stream
extern ProofStream* global_proof_stream;

uint64_t ***commit_host(Fri *fri, uint64_t **codeword, int codeword_len, uint64_t **tree_layer) {
    cout << "=== Commit Phase (Host) ===" << endl;
    uint64_t ***codewords = (uint64_t ***)malloc(fri_num_rounds(fri) + 1 * sizeof(uint64_t **));
    //get random bytes from Fiat-Shamir proof stash
    size_t num_bytes_alpha = 24; //assuming alpha fits in 24 bytes
    unsigned char* alpha_bytes = prover_fiat_shamir(global_proof_stream, num_bytes_alpha);
    //use libiop to sample field element from bytes
    //libiop::field_subset<FieldT> can sample random elements
    FieldT alpha = FieldT::random_element();
    //serialize alpha to bytes
    unsigned char* alpha_bytes = serialize_field_element(alpha);
    push_object(global_proof_stream, (void*)alpha_bytes);
    //allocate space for root hash (will be computed on GPU)
    memset(root_hash, 0, HASH_WORDS * 8);
    memset(alpha_hash_output, 0, 32);
    cout << "Alpha sampled and root hash space allocated." << endl;
    cout << "Codeword computations will happen on GPU kernel." << endl;
    for(int r = 0; r < fri->num_rounds; r++) {
        cout << "Round " << r << " processing..." << endl;
        //hash alpha for every round 
        hash_sha3_256(alpha, 24, alpha_hash_output); //assuming alpha fits in 24 bytes
        memcpy(alpha, alpha_hash_output, 24); //update alpha for next round
        cout << "Alpha hashed for round " << r << "." << endl;
        for(int i = 0; i < HASH_WORDS; i++) {
        cout << "Alpha hash output word " << i << ": " << ((uint64_t*)alpha_hash_output)[i] << endl;
        }
        bool is_last_round = (r == fri->num_rounds - 1);
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
        commit_launch(codeword, next_codeword, alpha, N, root, tree_layer, tree_layer_nxt, tree, last_round, is_last_round)
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
    push_object(global_proof_stream, (void*)root_hash);
    push_count++;
    cout << "  Pushed root hash" << endl;
    
    // 2. Push alpha
    push_object(global_proof_stream, (void*)alpha_bytes);
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
    cout << "  Pushed " << final_codeword_length << " last codeword elements" << endl;
    
    cout << "Total push count in commit: " << push_count << endl;
    
    free(alpha_bytes);
    return codewords;
}

size_t* prove(Fri* fri, std::vector<FieldT>& codeword, std::vector<std::vector<FieldT>>& merkle_tree_layers) {
    size_t codeword_length = codeword.size();
    cout << "codeword length in prove(): " << codeword_length << endl;
    
    // Create proof stream and run commit phase
    global_proof_stream = create_proof_stream();
    uint64_t ***codewords = commit_host(fri, codeword, codeword_length, merkle_tree_layers);
    
    // Sample indices using Fiat-Shamir (deterministic from transcript)
    size_t num_bytes = 32; // SHA3-256 hash size
    uint8_t *seed = prover_fiat_shamir(global_proof_stream, num_bytes);
    
    // Allocate space for indices
    size_t *top_level_indices = new size_t[fri->num_colinearity_tests];
    size_t *reduced_indices = new size_t[fri->num_colinearity_tests];
    
    // Sample indices deterministically from seed
    // This matches the CUDA version's sample_indices function
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
    
    // Run query phase for each FRI round
    for (size_t r = 0; r < fri_num_rounds(fri) - 1; r++) {
        // Query this round with current indices
        indices = query(fri, codewords, 
                       codewords[r], codeword_length >> r,  // current codeword
                       codewords[r+1], codeword_length >> (r+1),  // next codeword
                       indices, r);
    }
    
    free(seed);
    free(reduced_indices);
    
    return top_level_indices;
}