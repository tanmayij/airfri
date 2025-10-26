#include "../include/verify.cuh"
#include "../include/fri_utils.cuh"
#include "../include/commit_kernel.cuh"
#include "../include/poly_eval.cuh"
#include "../include/merkle.cuh"
#include "../include/fs_transform.cuh"
#include "../include/prove.cuh"
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cassert>
#include <cstring>

int verify(Fri *fri, uint64_t **polynomial_values, int degree) { 
    uint64_t one[FIELD_WORDS];
    field_one(one, FIELD_WORDS);
    const int max_fri_domains = sizeof(preon.fri_domains) / sizeof(preon.fri_domains[0]); //should be 16 for 256-bit params 
    int basis_len = fri_log_domain_length(fri); //basis length of initial codeword 
    uint64_t *eval_basis = populate_eval_basis(basis_len);
    if (eval_basis != NULL) {
        printf("Eval basis in commit phase for basis_len %d:\n", basis_len);
        for (int i = 0; i < basis_len; i++) {
            printf("%016lx ", eval_basis[i]);
        }
        printf("\n");
    }
    assert(eval_basis != NULL && "eval_basis are not correct");
    uint64_t offset = get_offset_for_basis(basis_len);

    if (offset != 0) {
        printf("Offset in commit phase for basis_len %d is: %016lx\n", basis_len, offset);
    } else {
        printf("Offset in commit phase not found or is zero.\n");
    } 
    assert(eval_basis != NULL && "eval_basis are not correct!");
    uint64_t *root_verify = (uint64_t *)malloc(HASH_SIZE * sizeof(uint64_t)); 
    size_t last_codeword_length = (int)pow(2, fri_log_domain_length(fri)) >> fri_num_rounds(fri); //should be 32 for FRI192

    //get the elements here too
    uint64_t var1[FIELD_WORDS] = {0x2}; // Represents var1
    uint64_t var2[FIELD_WORDS] = {0x3}; // Represents var2
    uint64_t elements[256*FIELD_WORDS];

    generate_elements(elements, var1, var2, field_words);
    uint64_t *alpha = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
    uint64_t *sampled_alpha = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
    uint8_t *alpha_bytes = verifier_fiat_shamir(global_proof_stream, 24);

    root_verify = (uint64_t *)pull_object(global_proof_stream);
    pull_count++;
    sampled_alpha = (uint64_t *)pull_object(global_proof_stream);
    pull_count++;

    printf("printing FRI Merkle root");
    for(int j = 0; j< FIELD_WORDS; j++){
        printf("%016lx ", root_verify[j]);
    }
    printf("\n");

    //get the last codeword elements 
    printf("last codeword length: %zu\n", last_codeword_length);
    uint64_t **last_codeword_arr = (uint64_t **)malloc(last_codeword_length * sizeof(uint64_t));
    for (int i = 0; i < last_codeword_length; i ++){
        *last_codeword_arr = (uint64_t *)malloc(field_words * sizeof(uint64_t *));
    }   
    //extract last codeword
    for(int i=0; i<last_codeword_length; i++){
        last_codeword_arr[i] = (uint64_t *)pull_object(global_proof_stream);
        pull_count++;
    }
    printf("\n");

    uint64_t *copied_last_codeword = (uint64_t *)malloc(HASH_SIZE);
    memcpy(copied_last_codeword, last_codeword_arr, HASH_SIZE);

    //degree fetch
    degree = (last_codeword_length * fri_num_rounds(fri) / fri->expansion_factor) - 1; //if last codeword length = 32, this should be 0
    //last basis and offset
    uint64_t *last_basis = (uint64_t *)malloc((int)log2(last_codeword_length) * sizeof(uint64_t));
    last_basis = populate_eval_basis((int)log2(last_codeword_length));
    uint64_t last_offset = get_offset_for_basis((int)log2(last_codeword_length));

    Domain last_domain;
    int r = fri_num_rounds(fri);

    initialize_domain(&last_domain, fri->initial_domain_length >> r, last_basis, &last_offset);
    uint64_t **last_domain_elements = (uint64_t **)malloc(last_domain.size * sizeof(uint64_t *));
    for (int i = 0; i < last_domain.size; i++) {
        last_domain_elements[i] = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t)); //this is the interpolant
        i_th_ele_in_span(last_domain_elements[i], last_basis, last_domain.basis_len, i);
        field_add(last_domain_elements[i], last_domain_elements[i], &last_offset, 1);
    }

    uint64_t poly_coeffs_[field_words * last_domain.size];
    interpolate_domain_single(poly_coeffs_, last_domain_elements, copied_last_codeword, last_domain.size, FIELD_WORDS);
    uint64_t eval_result[FIELD_WORDS];
    poly_eval_over_domain(eval_result, poly_coeffs_, last_domain.size, &last_domain, field_words, preon.field_bytesize); //this function calls fft
    if (memcmp(copied_last_codeword, eval_result, FIELD_WORDS * sizeof(uint64_t)) == 0) {
        printf("Re-evaluated codeword does not match original!\n");
        free(root_verify);
        for (int j = 0; j < last_domain.size; j++) {
            free(last_domain_elements[j]);
        }
        free(last_domain_elements);
        return 0;
    }
    if (poly_deg(poly_coeffs_, last_domain.size, FIELD_WORDS) <= degree) {
        printf("Low degree test passed\n");
        printf("Observed degree: %d\n", poly_deg(poly_coeffs_, last_domain.size, FIELD_WORDS) );
        printf("Should be no more than: %d\n", degree);
    } else {
        printf("Last codeword does not correspond to polynomial of low enough degree\n");
        printf("Observed degree: %d\n", poly_deg(poly_coeffs_, last_domain.size, FIELD_WORDS));
        printf("But should be no more than: %d\n", degree);
        return 0;
    } 
    for (r = 0; r < fri_num_rounds(fri) - 1; r++) {
        const int max_fri_domains = sizeof(preon.fri_domains) / sizeof(preon.fri_domains[0]);
        int basis_len = fri_log_domain_length(fri) - (r); 
        eval_basis = populate_eval_basis(basis_len);
        offset = get_offset_for_basis(basis_len);
        assert(eval_basis != NULL && "eval_basis are not correct!");

        //hash alpha
        hash_sha3_256((uint64_t *) sampled_alpha, 24, (uint64_t *) alpha); //hash alpha
        memcpy(sampled_alpha, alpha, 24); //repeat hash
        printf("Alpha that is hashed in verify: ");
        for (uint64_t  i = 0; i < field_words; i++) {
            printf("%016lx ", alpha[i]);
        }
        printf("\n");

        //get the indices
        size_t *c_indices = global_indices_tracker.round_indices[r];
        calculate_indices(c_indices, 
            global_query_indices.a_indices, 
            global_query_indices.b_indices, 
            fri->num_colinearity_tests);
        
        printf("Indices in verify:\n");
        for (int i = 0; i < fri->num_colinearity_tests; i++) {
            printf("a: %zu, b: %zu, c: %zu\n", global_query_indices.a_indices[i], global_query_indices.b_indices[i], c_indices[i]);
        }

        //arrays store information about codeword elements for each colinearity test
        uint64_t **aa = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **bb = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **cc = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **ay = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **by = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));
        uint64_t **cy = (uint64_t **)malloc(fri->num_colinearity_tests * sizeof(uint64_t *));

        // uint64_t precomputed_inverses[max_fri_domains] = {0}; //should i make this a public variable?
        // load_precomputed_inverses("fri_inverses_256.txt", precomputed_inverses);
        int exponent = preon.fri_domains[0].basis_len - fri_log_domain_length(fri); 
        uint64_t denominator_inv[field_words] = {0};
        
        denominator_inv[0] = precomputed_inverses[exponent + r]; 
        printf("denominator inverse %016lx\n", denominator_inv[0]);
        denominator_inv[1] = {0};
        denominator_inv[2] = {0};
        denominator_inv[3] = {0};
        //colinearity check. to change: ideally, we must check if merkle root verifies and then perform check. but for now this order is okay.
        for(int s = 0; s < fri -> num_colinearity_tests; s++){
            ay[s] = (uint64_t *)pull_object(global_proof_stream);
            pull_count++;
            print_field("ay", ay[s], field_words);
            aa[s] = ay[s];
            by[s] = (uint64_t *)pull_object(global_proof_stream);
            pull_count++;
            print_field("by", by[s], field_words);
            bb[s] = by[s];
            cy[s] = (uint64_t *)pull_object(global_proof_stream);
            pull_count++;
            print_field("cy", cy[s], field_words);
            cc[s] = cy[s];
         
            if (r == 0) {
                polynomial_values[global_query_indices.a_indices[s]] = aa[s];
                polynomial_values[global_query_indices.b_indices[s]] = bb[s];
            }
            // Printing the polynomial values that are not zero
            // printf("Polynomial values:\n");
            // int k = 0; // Define k outside the loop
            // for (int i = 0; i < fri->num_colinearity_tests; i++) {
            //     //check if the polynomial value is zero, skip if true
            //     if (is_zero(polynomial_values[i], 1)) {
            //         continue;
            //     }
            
            //     // Print the current non-zero polynomial value
            //     for (int j = 0; j < FIELD_WORDS; j++) {  //assuming FIELD_WORDS is the length of each polynomial value
            //         printf("%016lx ", polynomial_values[i][j]);
            //     }
            //     printf("\n");
            //     for (int j = 0; j < FIELD_WORDS; j++) {
            //         polynomial_values[k][j] = polynomial_values[i][j];
            //     }
            //     k++;
            // }
            uint64_t ax[field_words] = {0};
            uint64_t ax_temp[field_words] = {0};
            uint64_t bx[field_words] = {0};
            uint64_t bx_temp[field_words] = {0};
            uint64_t alpha_offset[field_words] = {0};
            uint64_t cx[field_words] = {0};
            uint64_t diff1[field_words] = {0};
            uint64_t diff2[field_words] = {0};
            uint64_t temp1[field_words] = {0};
            uint64_t temp2[field_words] = {0};
            uint64_t temp3[field_words] = {0};
            uint64_t temp4[field_words] = {0};
            uint64_t temp5[field_words] = {0};
            printf("basis len %d\n", basis_len);
            i_th_ele_in_span(ax_temp, eval_basis, basis_len, global_query_indices.a_indices[s]);
            field_add(ax, ax_temp, &offset, field_words);
            print_field("ax", ax, field_words);
            ax[1] = {0};
            ax[2] = {0};
            ax[3] = {0};
            i_th_ele_in_span(bx_temp, eval_basis, basis_len, global_query_indices.b_indices[s]);
            field_add(bx, &offset, bx_temp, field_words);
            print_field("bx", bx, field_words);
            bx[1] = {0};
            bx[2] = {0};
            bx[3] = {0};
            memcpy(cx, alpha, FIELD_WORDS * sizeof(uint64_t));
            field_sub(diff1, ay[s], by[s], field_words);
            print_field("diff1", diff1, field_words);
            field_sub(diff2, ax, bx, field_words);
            field_inv(temp1, diff2, field_words);
            memcpy(alpha_offset, ax, sizeof(uint64_t));
            field_sub(temp2, cx, ax, field_words);
            print_field("cx", cx, field_words);
            print_field("temp2", temp2, field_words);
            field_mul(temp3, temp2, denominator_inv, field_words);
            print_field("temp3", temp3, field_words);
            field_mul(temp4, temp3, diff1, field_words);
            print_field("temp4", temp4, field_words);
            field_add(temp5, temp4, ay[s], field_words);
            print_field("temp5", temp5, field_words);
            if (memcmp(cy[s], temp5, field_words * sizeof(uint64_t))!= 0) {
                printf("Colinearity check failed\n");
                printf("%d\n", memcmp(cy[s], temp5, field_words * sizeof(uint64_t)));
                free(aa);
                free(bb);
                free(cc);
                //free(root_verify);
                //free(alpha);
                for (int j = 0; j < last_domain.size; j++) {
                    free(last_domain_elements[j]);
                }
                free(last_domain_elements);
                return 0;
            }
            printf("Colinearity check passed");
            //check auth_path_a
            uint64_t **auth_path_a = (uint64_t **)malloc(MAX_PROOF_PATH_LENGTH * sizeof(uint64_t *));
            for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
                auth_path_a[i] = (uint64_t *)malloc(CONCAT_WORDS * sizeof(uint64_t));  // Allocate space for each hash
            }

            size_t *proof_len_ptr_a = (size_t *)pull_object(global_proof_stream);
            size_t proof_len_a = *proof_len_ptr_a;
            printf("proof len a here is: %zu \n", proof_len_a);
            for (size_t i = 0; i < proof_len_a; i++) {
                auth_path_a[i] = (uint64_t *)pull_object(global_proof_stream);
            }

            printf("Authentication Path:\n");
            for (size_t i = s; i < proof_len_a; i++) {
                printf("Layer %zu: ", i);

                size_t element_size;
                if (i == 0) {
                    element_size = FIELD_WORDS; //Layer 0: Only the sibling (FIELD_WORDS size)
                } else if (i < 13) {
                    element_size = FIELD_WORDS + FIELD_WORDS + HASH_WORDS; //Layers 1-12: Sibling includes FIELD_WORDS + HASH_WORDS
                } else {
                    element_size = HASH_WORDS; //Layers 13-16: Only the sibling hash (HASH_WORDS size)
                }

                for (size_t j = 0; j < element_size; j++) {
                    printf("%016lx ", auth_path_a[i][j]);
                }
                printf("\n");
            }
            printf("End of Authentication Path.\n");
            printf("pull in verify for auth_path_a %d\n", pull_count);
            if(r == 0) {
            if (!merkle_verify(root_verify, global_query_indices.a_indices[s], auth_path_a, proof_len_a, aa[s], r, r)) {
                printf("merkle authentication path verification fails for aa\n");
                free(aa);
                free(bb);
                free(cc);
                free(root_verify);
                for (int j = 0; j < last_domain.size; j++) {
                    free(last_domain_elements[j]);
                }
                free(last_domain_elements);
                return 0;
            }
            }

            // Free auth_path_a after use
            for (size_t i = 0; i < proof_len_a; i++) {
                free(auth_path_a[i]);
            }
            free(auth_path_a);
            //verify auth_path_b
            uint64_t **auth_path_b = (uint64_t **)malloc(MAX_PROOF_PATH_LENGTH * sizeof(uint64_t *));
            for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
                auth_path_b[i] = (uint64_t *)malloc(HASH_SIZE);  // Allocate space for each hash
            }
            size_t *proof_len_ptr_b = (size_t *)pull_object(global_proof_stream);
            size_t proof_len_b = *proof_len_ptr_b;
  
            for (size_t i = 0; i < proof_len_b; i++) {
                auth_path_b[i] = (uint64_t *)pull_object(global_proof_stream); // Pull each hash
            }
            pull_count++;
            printf("pull in verify for auth_path_b %d\n", pull_count);
            if(r == 0) {
            if (!merkle_verify(root_verify, global_query_indices.b_indices[s], auth_path_b, proof_len_b, bb[s], r, r)) {
                printf("merkle authentication path verification fails for bb\n");
                free(aa);
                free(bb);
                free(cc);
                //free(roots);
                //free(alphas);
                for (int j = 0; j < last_domain.size; j++) {
                    free(last_domain_elements[j]);
                }
                free(last_domain_elements);
                return 0;
            }
            }
            //verify auth_path_c
            uint64_t **auth_path_c = (uint64_t **)malloc(MAX_PROOF_PATH_LENGTH * sizeof(uint64_t *));
            for (size_t i = 0; i < MAX_PROOF_PATH_LENGTH; i++) {
                auth_path_c[i] = (uint64_t *)malloc(HASH_SIZE);  // Allocate space for each hash
            }
            size_t *proof_len_ptr_c = (size_t *)pull_object(global_proof_stream);
            size_t proof_len_c = *proof_len_ptr_c;

            for (size_t i = 0; i < proof_len_c; i++) {
                auth_path_c[i] = (uint64_t *)pull_object(global_proof_stream); // Pull each hash
            }
            pull_count++;
            printf("pull in verify for auth_path_b %d\n", pull_count);
            // if (!merkle_verify(root_verify, global_query_indices.c_indices[s], auth_path_c, proof_len_c, cc[s], r+1, r)) {
            //     printf("merkle authentication path verification fails for cc\n");
            //     free(aa);
            //     free(bb);
            //     free(cc);
            //     //free(roots);
            //     //free(alphas);
            //     for (int j = 0; j < last_domain.size; j++) {
            //         free(last_domain_elements[j]);
            //     }
            //     free(last_domain_elements);
            //     return 0;
            // }
        }
        free(global_query_indices.c_indices);
        free(global_query_indices.a_indices);
        free(global_query_indices.b_indices);
        free(aa);
        free(bb);
        free(cc);
        //free(precomputed_inverses);

    }

    printf("Verification successful!\n");
    return 1;

}