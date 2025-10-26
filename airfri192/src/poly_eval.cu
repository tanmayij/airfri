#include "../include/poly.cuh"

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <stdbool.h>

#include "../include/domain.cuh"
//#include "finite_fieldelement.h"
#include "../include/field.cuh"
//#include "fft.h"
//#include "fri.cuh"
// #include "additivefft.h"
// Polynomial structure
typedef struct {
    uint64_t **coefficients; // Array of field elements
    int degree;              // Degree of the polynomial
} Polynomial;


int poly_add(uint64_t *sum, const uint64_t *a, const size_t a_len, const uint64_t *b, const size_t b_len, const size_t field_words)
{
    size_t sum_len = a_len;
    if (b_len > a_len)
        sum_len = b_len;
    // Assume sum is already allocated?
    // sum = realloc(ptr, new_size);
    // new_ptr = realloc(ptr, new_size);
    // if(new_ptr == null){  //用 !new_ptr 檢查也可以
    //     // 錯誤處理 error handling
    // }
    // ptr = new_ptr
    uint64_t *zero = (uint64_t *)malloc(field_words * sizeof(uint64_t));
    memset(zero, 0, field_words * sizeof(uint64_t));
    for (size_t i = 0; i < sum_len; i++)
    {
        if (i < a_len && i < b_len)
            field_add(&sum[i * field_words], &a[i * field_words], &b[i * field_words], field_words);
        else if (i >= a_len && i < b_len)
            field_add(&sum[i * field_words], zero, &b[i * field_words], field_words);
        else if (i >= b_len && i < a_len)
            field_add(&sum[i * field_words], &a[i * field_words], zero, field_words);
        else
            return -1; // Or any error code about out of range
    }
    free(zero);
    return 0;
}

void poly_addEqual(uint64_t *a, const size_t a_len, const uint64_t *b, const size_t b_len, const size_t field_words)
{
    // A potential reallocation could happen here
    // Or is it ok to only support a_len >= b_len?
    // So that a is already allocated to desired length
    assert(a_len >= b_len);
    for (size_t i = 0; i < a_len * field_words; i += field_words)
    {
        if (i < b_len * field_words)
            field_addEqual(&a[i], &b[i], field_words);
    }
}

void poly_scalarMul(uint64_t *product, const uint64_t *a, const size_t a_len, const uint64_t *b, const size_t field_words)
{
    memset(product, 0, a_len * field_words * sizeof(uint64_t));
    if (is_zero(b, field_words))
    {
        return;
    }

    for (size_t i = 0; i < a_len * field_words; i += field_words)
    {
        if (!is_zero(&a[i], field_words))
        {
            field_mul(&product[i], &a[i], b, field_words);
        }
    }
}

void poly_scalarMulEqual(uint64_t *a, const size_t a_len, const uint64_t *b, const size_t field_words)
{
    if (is_zero(b, field_words))
    {
        memset(a, 0, a_len * field_words * sizeof(uint64_t));
        return;
    }
    for (size_t i = 0; i < a_len * field_words; i += field_words)
    {
        if (!is_zero(&a[i], field_words))
        {
            field_mulEqual(&a[i], b, field_words);
        }
    }
}

void poly_mul(uint64_t *product, const uint64_t *a, const size_t a_len, const uint64_t *b, const size_t b_len, const size_t field_words)
{
    // product should be as length a_len + b_len - 1
    size_t product_len = a_len + b_len - 1;
    memset(product, 0, product_len * field_words * sizeof(uint64_t));
    uint64_t *temp = (uint64_t *)malloc(b_len * field_words * sizeof(uint64_t));
    memset(temp, 0, b_len * field_words * sizeof(uint64_t));
    for (size_t i = 0; i < a_len * field_words; i += field_words)
    {
        if (!is_zero(&a[i], field_words))
        {
            poly_scalarMul(temp, b, b_len, &a[i], field_words);
            poly_addEqual(&product[i], b_len, temp, b_len, field_words);
        }
    }
    free(temp);
}
void poly_mulEqual(uint64_t *a, const size_t a_len, const uint64_t *b, const size_t b_len, const size_t field_words)
{
    // May not be as straightforward as poly_mul
    // Huge cost when duplicating, probably not implementing this
    assert("not implemented yet" && 0);
}

void poly_div(uint64_t *quotient, uint64_t *remainder, const uint64_t *a, const size_t a_len, const uint64_t *b, const size_t b_len, const size_t field_words)
{
    // a = b * quotient + remainder
    // inverse of the leading term of b
    size_t field_bytesize = field_words * sizeof(uint64_t);
    uint64_t *leadingTermInverse = (uint64_t *)malloc(field_bytesize);
    field_inv(leadingTermInverse, &b[(b_len - 1) * field_words], field_words);
    const size_t remainder_len = b_len - 1;

    if (a_len < b_len)
    {
        memset(quotient, 0, field_bytesize);
        memcpy(remainder, a, a_len * field_bytesize);
        free(leadingTermInverse);
        return;
    }

    memcpy(remainder, a, remainder_len * field_bytesize);
    memcpy(quotient, &a[remainder_len * field_words], (a_len - remainder_len) * field_bytesize);

    for (size_t i = (a_len - remainder_len) * field_words; i > 0; i -= field_words)
    {
        // printf("%zu\n", i);
        uint64_t *twist = (uint64_t *)malloc(field_bytesize);
        field_mul(twist, &quotient[i - field_words], leadingTermInverse, field_words);
        memcpy(&quotient[i - field_words], twist, field_bytesize);

        /* subtract twist*Z * y^i thus clearing the i-th term of P */
        if (b_len >= 2 && !is_zero(twist, field_words))
        {
            for (size_t j = (b_len - 1) * field_words; j > 0; j -= field_words)
            {
                // printf("%zu\n", j);
                uint64_t *deduction = (uint64_t *)malloc(field_bytesize);
                field_mul(deduction, twist, &b[j - field_words], field_words);
                if (!is_zero(deduction, field_words))
                {
                    if ((i - field_words + j - field_words) < remainder_len * field_words)
                        field_addEqual(&remainder[i - field_words + j - field_words], deduction, field_words);
                    else
                        field_addEqual(&quotient[i - field_words + j - field_words - remainder_len * field_words], deduction, field_words);
                }
                free(deduction);
            }
        }
        free(twist);
    }
    free(leadingTermInverse);
}


void poly_eval(uint64_t *result, const uint64_t *poly, const size_t poly_len, const uint64_t *x, const size_t field_words)
{
    memset(result, 0, field_words * sizeof(uint64_t));
    uint64_t term[field_words];
    for (size_t i = 0; i < poly_len; i++)
    {
        if (i == 0)
        {
            memcpy(term, &poly[i * field_words], field_words * sizeof(uint64_t));
        }
        else
        {
            field_mul(term, &poly[i * field_words], x, field_words);
            for (size_t j = 1; j < i; j++)
            {
                field_mulEqual(term, x, field_words);
            }
        }
        field_addEqual(result, term, field_words);
    }
}

void poly_eval_over_domain(uint64_t *result, const uint64_t *poly, const size_t poly_len, const Domain *domain, const size_t field_words, const size_t field_bytesize)
{
    fft(result, poly, poly_len, field_words, domain->size, domain->shift);
    // TODO: fix below optimized code
    /*
        We implement /affine/ linearized polynomials for which the
        bilinear property is not directly applicable, as we need to make
        sure that constant term is included only once.

        Therefore, evaluating over subspace below, we subtract constant
        term from evaluations over the basis, but include the constant
        term in the shift calculation.
        */

    // size_t domain_dimension = domain->basis_len;
    // uint64_t *eval_at_basis = (uint64_t *)malloc(domain_dimension * field_bytesize);
    // for (size_t i = 0; i < domain_dimension; i++)
    // {
    //     poly_eval(&eval_at_basis[i * field_words], poly, poly_len, &domain->basis[i * field_words], field_words);
    //     field_subEqual(&eval_at_basis[i * field_words], &poly[0], field_words);
    // }
    // uint64_t *eval_at_shift = (uint64_t *)malloc(field_bytesize);
    // poly_eval(eval_at_shift, poly, poly_len, domain->shift, field_words);
    // all_subset_sums(result, eval_at_basis, domain_dimension, eval_at_shift, field_words);

    // free(eval_at_basis);
    // free(eval_at_shift);
}

bool is_zero_poly(const uint64_t *coefficients, size_t len, size_t field_words) {
    uint64_t *zero = (uint64_t *)malloc(field_words * sizeof(uint64_t));
    memset(zero, 0, field_words * sizeof(uint64_t)); // Create a zero element

    for (size_t i = 0; i < len; i++) {
        if (memcmp(&coefficients[i * field_words], zero, field_words * sizeof(uint64_t)) != 0) {
            free(zero);
            return false; // Found a non-zero coefficient
        }
    }

    free(zero);
    return true; // All coefficients are zero
}

int poly_deg(const uint64_t *poly, size_t poly_len, size_t field_words) {
    if (is_zero_poly(poly, poly_len, field_words)) {
        return -1; // Polynomial is zero
    }

     size_t maxindex = 0;
    uint64_t *zero = (uint64_t *)malloc(field_words * sizeof(uint64_t));
    memset(zero, 0, field_words * sizeof(uint64_t)); // Create a zero element

    for (size_t i = 0; i < poly_len; i++) {
        if (memcmp(&poly[i * field_words], zero, field_words * sizeof(uint64_t)) != 0) {
            maxindex = i * field_words; // Update maxindex to the current index if the coefficient is non-zero
        }
    }

    free(zero);
    return maxindex / (field_words * sizeof(uint64_t)); 
}

void print_field(const char *label, const uint64_t *field, size_t field_words);

void interpolate_domain(uint64_t *result, uint64_t **domain, uint64_t **values, size_t domain_length, size_t field_words)
{
    assert(domain_length > 0 && "cannot interpolate between zero points");
    memset(result, 0, domain_length * field_words * sizeof(uint64_t));
    uint64_t *temp_poly = (uint64_t *)malloc(field_words * sizeof(uint64_t));
    uint64_t *prod = (uint64_t *)malloc(field_words * sizeof(uint64_t));
    uint64_t *diff = (uint64_t *)malloc(field_words * sizeof(uint64_t));
    uint64_t *inv_diff = (uint64_t *)malloc(field_words * sizeof(uint64_t));
    uint64_t *temp_term = (uint64_t *)malloc(field_words * sizeof(uint64_t));

    for (size_t i = 0; i < domain_length; i++)
    {
        memcpy(temp_poly, values[i], field_words * sizeof(uint64_t));

        for (size_t j = 0; j < domain_length; j++)
        {
            if (i != j)
            {
                // Calculate the inverse difference
                field_sub(diff, domain[i], domain[j], field_words);
                //printf("domain[i]: %016llx domain[j]: %016llx\n", domain[i], domain[j]);
                field_inv(inv_diff, diff, field_words);

                // Initialize prod to one
                memset(prod, 0, field_words * sizeof(uint64_t));
                prod[field_words - 1] = 1;

                for (size_t k = 0; k < domain_length; k++)
                {
                    if (k != j)
                    {
                        // Calculate term (x - domain[k])
                        memset(temp_term, 0, field_words * sizeof(uint64_t));
                        field_neg(temp_term, domain[k], field_words);

                        // Multiply term by inv_diff and accumulate in prod
                        poly_scalarMulEqual(temp_term, domain_length, inv_diff, field_words);
                        poly_mul(prod, prod, domain_length, temp_term, domain_length, field_words);
                    }
                }
                // Multiply temp_poly by the accumulated product
                poly_scalarMulEqual(temp_poly, domain_length, prod, field_words);
            }
        }
        // Accumulate the result
        poly_addEqual(result, domain_length, temp_poly, 1, field_words);
    }

    free(temp_poly);
    free(prod);
    free(diff);
    free(temp_term);
}


void interpolate_domain_single(uint64_t *result, uint64_t **domain, uint64_t *values, size_t domain_length, size_t field_words) {
    // Assertions
    assert(domain_length > 0 && "cannot interpolate between zero points");

    // Initialize necessary variables
    size_t x_len = 2 * field_words; // x will contain 2 field elements
    uint64_t *x = (uint64_t *)malloc(x_len * sizeof(uint64_t));

    memset(x, 0, x_len * sizeof(uint64_t));
    x[(2 * field_words) - 1] = 1; // Set x to [0, 1]

    // for (size_t k = 0; k < 2 * field_words; k++) {
    //     printf("%016llx ", x[k]);
    // }

    uint64_t *acc = (uint64_t *)malloc(domain_length * field_words * sizeof(uint64_t));
    memset(acc, 0, domain_length * field_words * sizeof(uint64_t));

    uint64_t *prod = (uint64_t *)malloc(domain_length * field_words * sizeof(uint64_t));
    uint64_t *prod_temp = (uint64_t *)malloc(domain_length * field_words * sizeof(uint64_t));
    uint64_t *term = (uint64_t *)malloc(field_words * sizeof(uint64_t));
    uint64_t *diff = (uint64_t *)malloc(field_words * sizeof(uint64_t));
    uint64_t *inv_diff = (uint64_t *)malloc(field_words * sizeof(uint64_t));

    //printf("Starting interpolation with domain length: %zu and field words: %zu\n", domain_length, field_words);

    for (size_t i = 0; i < domain_length; i++) {
        // prod = Polynomial([values[i]])
        memcpy(prod, &values[i * field_words], field_words * sizeof(uint64_t));
        // printf("Initial prod for i=%zu: ", i);
        // for (size_t k = 0; k < field_words; k++) {
        //     printf("%016llx ", prod[k]);
        // }
        // printf("\n");

        for (size_t j = 0; j < domain_length; j++) {
            if (j == i) continue;

            // term = x - Polynomial([domain[j]])
            memcpy(term, x, x_len * sizeof(uint64_t));
            // printf("domain [0] : ");
            // for (size_t k = 0; k < 1; k++) {
            //     printf("%016llx ", *domain[0]);
            // }
            field_sub(&term[0], &term[0], domain[j], field_words);
    

            // prod = prod * (x - Polynomial([domain[j]]))
            memset(prod_temp, 1, field_words * sizeof(uint64_t));
            poly_mul(prod_temp, prod, 1, term, 1, field_words);
            // printf("prod after multiplication for i=%zu, j=%zu: ", i, j);
            // for (size_t k = 0; k < field_words; k++) {
            //     printf("%016llx ", prod_temp[k]);
            // }
            // printf("\n");

            // Calculate the inverse of (domain[i] - domain[j])
            field_sub(diff, domain[i], domain[j], field_words);
            field_inv(inv_diff, diff, field_words);
            // printf("diff and inv_diff for i=%zu, j=%zu: ", i, j);
            // for (size_t k = 0; k < field_words; k++) {
            //     printf("diff=%016llx inv_diff=%016llx ", diff[k], inv_diff[k]);
            // }
            // printf("\n");

            // Multiply prod by the inverse difference
            poly_scalarMulEqual(prod_temp, domain_length, inv_diff, field_words);
            // printf("prod after scalar multiplication for i=%zu, j=%zu: ", i, j);
            // for (size_t k = 0; k < field_words; k++) {
            //     printf("%016llx ", prod[k]);
            // }
            // printf("\n");
        }
        // acc = acc + prod
        poly_addEqual(acc, domain_length, prod_temp, 1, field_words);
        // printf("acc after addition for i=%zu: ", i);
        // for (size_t k = 0; k < domain_length * field_words; k++) {
        //     printf("%016llx ", acc[k]);
        // }
        // printf("\n");
    }

    // Copy result
    memcpy(result, acc, field_words * sizeof(uint64_t));
    // printf("Final result: ");
    // for (size_t k = 0; k < field_words; k++) {
    //     printf("%016llx ", result[k]);
    // }
    // printf("\n");

    // Free allocated memory - this throws an error, figure out why.
    // free(x);
    // free(acc);
    // free(prod);
    // free(term);
    // free(diff);
    // free(inv_diff);
}
// bool test_colinearity(uint64_t **points, size_t num_points, size_t field_words)
// {
//     uint64_t **domain = (uint64_t **)malloc(num_points * sizeof(uint64_t *));
//     uint64_t **values = (uint64_t **)malloc(num_points * sizeof(uint64_t *));
//     for (size_t i = 0; i < num_points; i++)
//     {
//         domain[i] = points[i * 2];
//         values[i] = points[i * 2 + 1];
//     }

//     uint64_t *interpolated_poly = (uint64_t *)malloc(num_points * field_words * sizeof(uint64_t));
//     interpolate_domain(interpolated_poly, domain, values, num_points, field_words);

//     size_t degree = poly_deg(interpolated_poly, num_points, field_words);

//     free(interpolated_poly);
//     free(domain);
//     free(values);

//     return degree == 1;
// }

void shift_poly(uint64_t *result, const uint64_t *poly, size_t poly_len, const uint64_t *offset, size_t field_words)
{
    uint64_t *xi = (uint64_t *)malloc(field_words * sizeof(uint64_t));
    uint64_t *temp_poly = (uint64_t *)malloc((poly_len + 1) * field_words * sizeof(uint64_t));
    memset(result, 0, (poly_len + 1) * field_words * sizeof(uint64_t));
    memset(xi, 0, field_words * sizeof(uint64_t));
    memcpy(xi, offset, field_words * sizeof(uint64_t));

    for (size_t i = 0; i < poly_len; i++)
    {
        memcpy(temp_poly, poly + i * field_words, field_words * sizeof(uint64_t));
        for (size_t j = 0; j < i; j++)
        {
            field_mulEqual(temp_poly, xi, field_words);
        }
        poly_addEqual(result, poly_len + 1, temp_poly, 1, field_words);
    }

    free(xi);
    free(temp_poly);
}

void scale_poly(uint64_t *result, const uint64_t *poly, size_t poly_len, const uint64_t *factor, size_t field_words)
{
    memset(result, 0, poly_len * field_words * sizeof(uint64_t));
    uint64_t *factor_power = (uint64_t *)malloc(field_words * sizeof(uint64_t));
    memcpy(factor_power, factor, field_words * sizeof(uint64_t));

    for (size_t i = 0; i < poly_len; i++)
    {
        memcpy(&result[i * field_words], &poly[i * field_words], field_words * sizeof(uint64_t));
        for (size_t j = 0; j < i; j++)
        {
            field_mulEqual(&result[i * field_words], factor_power, field_words);
        }
    }

    free(factor_power);
}
// #include <cuda_runtime.h>
// #include <stdint.h>
// #include <stdio.h>
// #include <cstring>
// #include "poly_eval.hpp"

// #define FIELD_WORDS 4

// // Forward declaration of field operations (you'll need to include or implement these)
// __device__ void field_mul(uint64_t *result, const uint64_t *a, const uint64_t *b, size_t field_words);
// __device__ void field_mulEqual(uint64_t *a, const uint64_t *b, size_t field_words);
// __device__ void field_addEqual(uint64_t *a, const uint64_t *b, size_t field_words);

// // Kernel to perform poly_eval in parallel on the GPU
// __global__ void poly_eval_kernel(
//     uint64_t *codeword,        // Flattened 1D array
//     uint64_t *poly, 
//     int poly_len, 
//     uint64_t *domain_elements, // Flattened 1D array
//     size_t field_words, 
//     size_t initial_domain_length
// ) {
//     int i = blockIdx.x * blockDim.x + threadIdx.x;

//     if (i < initial_domain_length) {
//         // Compute the starting index for the i-th element in the flattened codeword array
//         uint64_t *result = &codeword[i * field_words];

//         // Initialize the result to zero
//         for (size_t k = 0; k < field_words; k++) {
//             result[k] = 0;
//         }

//         uint64_t term[FIELD_WORDS];  // Temporary array for term calculation

//         // Perform polynomial evaluation
//         for (size_t j = 0; j < poly_len; j++) {
//             if (j == 0) {
//                 // Copy poly[j] to term
//                 for (size_t k = 0; k < field_words; k++) {
//                     term[k] = poly[j * field_words + k];
//                 }
//             } else {
//                 // Compute the starting index for the i-th domain element
//                 uint64_t *domain_elem = &domain_elements[i * field_words];

//                 // Multiply term by poly[j] and domain_elements[i]
//                 field_mul(term, &poly[j * field_words], domain_elem, field_words);

//                 // Repeatedly multiply by domain_elements[i] if necessary
//                 for (size_t k = 1; k < j; k++) {
//                     field_mulEqual(term, domain_elem, field_words);
//                 }
//             }

//             // Accumulate the result
//             field_addEqual(result, term, field_words);
//         }
//     }
// }

// // CUDA wrapper function - C linkage for calling from C++
// extern "C" void parallel_poly_eval_cuda(
//     uint64_t **codeword, 
//     uint64_t *poly_coeffs, 
//     int poly_len, 
//     uint64_t **domain_elements, 
//     size_t field_words, 
//     size_t initial_domain_length
// ) {
//     // GPU memory allocations
//     uint64_t *d_codeword, *d_poly_coeffs, *d_domain_elements;

//     // Flatten the codeword and domain_elements on the host
//     uint64_t *flattened_codeword = (uint64_t *)malloc(initial_domain_length * field_words * sizeof(uint64_t));
//     uint64_t *flattened_domain_elements = (uint64_t *)malloc(initial_domain_length * field_words * sizeof(uint64_t));

//     // Flatten codeword
//     for (size_t i = 0; i < initial_domain_length; ++i) {
//         for (size_t j = 0; j < field_words; ++j) {
//             flattened_codeword[i * field_words + j] = codeword[i][j];
//         }
//     }

//     // Flatten domain_elements
//     for (size_t i = 0; i < initial_domain_length; ++i) {
//         for (size_t j = 0; j < field_words; ++j) {
//             flattened_domain_elements[i * field_words + j] = domain_elements[i][j];
//         }
//     }

//     // Allocate memory on the GPU for codeword, poly_coeffs, and domain_elements
//     cudaMalloc((void **)&d_codeword, initial_domain_length * field_words * sizeof(uint64_t));
//     cudaMalloc((void **)&d_poly_coeffs, poly_len * field_words * sizeof(uint64_t));
//     cudaMalloc((void **)&d_domain_elements, initial_domain_length * field_words * sizeof(uint64_t));

//     // Copy data from host to device
//     cudaMemcpy(d_codeword, flattened_codeword, initial_domain_length * field_words * sizeof(uint64_t), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_poly_coeffs, poly_coeffs, poly_len * field_words * sizeof(uint64_t), cudaMemcpyHostToDevice);
//     cudaMemcpy(d_domain_elements, flattened_domain_elements, initial_domain_length * field_words * sizeof(uint64_t), cudaMemcpyHostToDevice);

//     // Define the number of threads and blocks
//     int threads_per_block = 256;
//     int num_blocks = ((initial_domain_length + threads_per_block - 1) / threads_per_block) * 2; //512 because 512 * 256 = 131072

//     // Launch the kernel
//     poly_eval_kernel<<<num_blocks, threads_per_block>>>(d_codeword, d_poly_coeffs, poly_len, d_domain_elements, field_words, initial_domain_length);

//     // Wait for the kernel to finish
//     cudaDeviceSynchronize();

//     // Copy the results back to the host (unflatten)
//     cudaMemcpy(flattened_codeword, d_codeword, initial_domain_length * field_words * sizeof(uint64_t), cudaMemcpyDeviceToHost);

//     // Unflatten the codeword back into the original 2D array
//     for (size_t i = 0; i < initial_domain_length; ++i) {
//         for (size_t j = 0; j < field_words; ++j) {
//             codeword[i][j] = flattened_codeword[i * field_words + j];
//         }
//     }

//     // Free GPU memory
//     cudaFree(d_codeword);
//     cudaFree(d_poly_coeffs);
//     cudaFree(d_domain_elements);

//     // Free host flattened arrays
//     free(flattened_codeword);
//     free(flattened_domain_elements);
// }

// // C++ wrapper that converts libff::gf256 vectors to raw arrays
// namespace poly_eval {

// std::vector<libff::gf256> evaluate_polynomial_on_domain_gpu(
//     const std::vector<libff::gf256>& poly_coeffs,
//     const std::vector<libff::gf256>& domain_elements
// ) {
//     size_t poly_len = poly_coeffs.size();
//     size_t domain_len = domain_elements.size();
//     size_t field_words = FIELD_WORDS;

//     // Allocate 2D arrays for codeword and domain elements
//     std::vector<std::vector<uint64_t>> codeword(domain_len, std::vector<uint64_t>(field_words, 0));
//     std::vector<std::vector<uint64_t>> domain_raw(domain_len, std::vector<uint64_t>(field_words, 0));
//     std::vector<uint64_t> poly_raw(poly_len * field_words, 0);

//     // Pack poly_coeffs into poly_raw
//     for (size_t i = 0; i < poly_len; ++i) {
//         const libff::gf256& elem = poly_coeffs[i];
//         elem.as_words(&poly_raw[i * field_words]);
//     }

//     // Pack domain_elements into domain_raw
//     for (size_t i = 0; i < domain_len; ++i) {
//         domain_elements[i].as_words(domain_raw[i].data());
//     }

//     // Prepare pointers for CUDA call
//     std::vector<uint64_t*> codeword_ptrs(domain_len);
//     std::vector<uint64_t*> domain_ptrs(domain_len);
//     for (size_t i = 0; i < domain_len; ++i) {
//         codeword_ptrs[i] = codeword[i].data();
//         domain_ptrs[i] = domain_raw[i].data();
//     }

//     // Call CUDA kernel
//     parallel_poly_eval_cuda(
//         codeword_ptrs.data(),
//         poly_raw.data(),
//         static_cast<int>(poly_len),
//         domain_ptrs.data(),
//         field_words,
//         domain_len
//     );

//     // Unpack codeword back to libff::gf256
//     std::vector<libff::gf256> result(domain_len);
//     for (size_t i = 0; i < domain_len; ++i) {
//         result[i] = libff::gf256::from_words(codeword[i].data());
//     }
//     return result;
// }

// } // namespace poly_eval
