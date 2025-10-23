#include <iostream>
#include <cmath>
#include <algorithm>
#include <ctime>
#include <vector>
#include <libff/algebra/fields/binary/gf256.hpp>
#include <libiop/algebra/field_subset/subspace.hpp>
#include <libiop/algebra/utils.hpp>
#include <libiop/algebra/fft.hpp>
#include "Cantor/cantor_basis.hpp"
#include "poly_eval.hpp"
#include "src/prove.hpp"

using namespace std;

//constants
constexpr size_t FIELD_WORDS = 4; //gf256 uses 4 uint64_t words
constexpr size_t HASH_WORDS = 4;  //sha3-256 hash is 32 bytes = 4 uint64_t words

void save_merkle_tree_to_file(uint64_t ***tree, int max_layers, int concat_words, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file for writing\n");
        return;
    }

    fprintf(file, "========= FULL MERKLE TREE DUMP =========\n");

    for (int layer = 0; layer < max_layers; layer++) {
        if(layer == 0 || (layer > 12)) { concat_words = 4; }
        else { concat_words = 8; }
        if (tree[layer] == NULL) {
            fprintf(file, "Layer %d is NULL\n", layer);
            continue;
        }

        int num_nodes = 1 << (max_layers - (layer)); // Approximate number of nodes

        fprintf(file, "\n=== Layer %d (%d elements) ===\n", layer, num_nodes);
        for (int node = 0; node < num_nodes; node++) {
            if (tree[layer][node] == NULL) {
                fprintf(file, "  Node %d: NULL\n", node);
                continue;
            }

            fprintf(file, "  Node %d: ", node);
            for (int j = 0; j < concat_words; j++) {
                fprintf(file, "%016lx ", tree[layer][node][j]);
            }
            fprintf(file, "\n");
        }
    }

    fprintf(file, "\n=========================================\n");

    fclose(file);
    printf("Merkle tree has been saved to %s\n", filename);
}

Fri* init_fri(int initial_domain_length, int expansion_factor, int num_colinearity_tests) {
    Fri* fri = new Fri;
    fri->initial_domain_length = initial_domain_length;
    fri->expansion_factor = expansion_factor;
    fri->num_colinearity_tests = num_colinearity_tests;
    return fri;
}

void test_fri(){
    //terms are hardcoded for GF(2^256) field and 192-bit security
    int degree = 32767; // 2^15 - 1
    int expansion_factor = 32;
    int num_colinearity_tests = 26;

    int initial_domain_length = (degree + 1) * expansion_factor; //2^20
    cout << "Initial Domain Length: " << initial_domain_length << std::endl;
    int num_rounds = static_cast<int>(log2(initial_domain_length));
    //use cantor basis
    typedef libff::gf256 FieldT;
    const size_t m = 20; //use first 20 basis elements
    std::vector<FieldT> basis = cantor_basis<FieldT>(m);
    
    //reverse the order of basis elements
    // std::vector<FieldT> basis(full_basis.begin(), full_basis.begin() + 20);
    // std::reverse(basis.begin(), basis.end());
    cout << "basis (first 20 rows):" << std::endl;
    for (size_t i = 0; i < basis.size(); ++i) {
        cout << "basis[" << i << "] = " << basis[i] << std::endl;
    }

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> result;
    assert(domain.size() != 0);
    Fri *fri = init_fri(initial_domain_length, expansion_factor, num_colinearity_tests);
    assert(fri != nullptr);
    int num_rounds = static_cast<int>(log2(initial_domain_length));
    cout << "number of rounds: " << num_rounds << endl;
    //generate random polynomial coefficients
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(degree + 1);
    //evaluate polynomial on domain to get codewords
    //use available hardware threads minus 1, or default to 31
    int num_threads = std::thread::hardware_concurrency() > 0 ? 
                  std::thread::hardware_concurrency() - 1 : 31;
    std::vector<FieldT> codeword = cantor::additive_FFT(poly_coeffs, m, num_threads);
    // std::vector<FieldT> codeword = poly_eval::evaluate_polynomial_on_domain(
    // poly_coeffs, domain.all_elements());

    cout << "codeword length: " << codeword.size() << endl;

    //initialize merkle tree with 21 layers (0-20)
    //layer 0: 2^20 leaves, layer 20: 1 root
    //FRI: 15 rounds (2^20 -> 2^5), so layers 0-15 have field elements
    //integrated structure: layers 0-15 have field elements, 16-20 hashes only
    std::vector<std::vector<FieldT>> merkle_tree_layers(21);
    //layer 0: store first codeword elements (2^20 elements)
    merkle_tree_layers[0] = codeword;
    
    cout << "initialized layer 0 with " << merkle_tree_layers[0].size() << " codeword elements" << endl;
    
    // Allocate storage for queried polynomial evaluation points
    // These will be filled during verify() with values at a_indices and b_indices
    // Size is 2 * initial_domain_length to accommodate indices from first folding round
    // In practice, only 2 * num_colinearity_tests (52) entries will be used
    const size_t points_size = initial_domain_length / 2; // 2^19 for safety
    std::vector<FieldT> points(points_size);
    
    cout << "allocated points array with " << points_size << " entries" << endl;
    cout << "testing valid fri instance" << endl;
    clock_t t_prover_start = clock();
    prove(fri, codeword, merkle_tree_layers);
    clock_t t_prover_end = clock();
    double prover_time = double(t_prover_end - t_prover_start) / CLOCKS_PER_SEC;
    cout << "Prover time: " << prover_time << " seconds" << endl;

    //save merkle tree to file for debugging
    save_merkle_tree_to_file(merkle_tree_layers, 21, 8, "../data/merkle_tree_dump.txt");

    // Run verifier
    clock_t t_verifier_start = clock();
    int verdict = verify(fri, points.data(), degree);
    clock_t t_verifier_end = clock();
    double verifier_time = double(t_verifier_end - t_verifier_start) / CLOCKS_PER_SEC;
    
    cout << "Verifier time: " << verifier_time << " seconds" << endl;

    if (verdict == 1) {
        cout << "Proof accepted" << endl;
    } else {
        cout << "Proof rejected (but should be valid!)" << endl;
        return;
    }
    
    // Final polynomial evaluation check
    // Verify that the polynomial from last codeword evaluates correctly on initial domain
    cout << "Performing final polynomial evaluation check..." << endl;
    
    // For each colinearity test, verify polynomial evaluation
    for (int i = 0; i < num_colinearity_tests; i++) {
        // Compute i-th element in span of Cantor basis
        // domain[i] = offset + sum(basis[j] where bit j of i is set)
        FieldT domain_point = FieldT::zero();
        
        for (int bit = 0; bit < m; bit++) {
            if (i & (1 << bit)) {
                domain_point += basis[bit];
            }
        }
        
        // Add offset (shift of affine subspace)
        domain_point += domain.shift();
        
        // Evaluate polynomial at this domain point
        // poly_eval: result = sum(poly_coeffs[j] * domain_point^j)
        FieldT eval_result = FieldT::zero();
        FieldT domain_power = FieldT::one();
        for (size_t j = 0; j < poly_coeffs.size(); j++) {
            eval_result += poly_coeffs[j] * domain_power;
            domain_power *= domain_point;
        }
        
        // Get the value from points array (filled during verify)
        // Points at even indices are x-coords, odd indices are y-coords
        FieldT expected_y = points[2 * i + 1];
        
        // Verify it matches
        if (eval_result != expected_y) {
            cout << "✗ Polynomial evaluates to wrong value at point " << i << endl;
            cout << "Expected: " << expected_y << endl;
            cout << "Got: " << eval_result << endl;
            assert(false);
        }
    }
    
    cout << "All polynomial evaluations correct!" << endl;
    cout << "Success! \\o/" << endl;
    

}
int main() {
    // libff::gf256::init_field();
    
    try {
        test_fri();
        cout << "test completed successfully" << endl;
    } catch (const std::exception& e) {
        cerr << "error: " << e.what() << endl;
        return 1;
    }
    
    return 0;
}
