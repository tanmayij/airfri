#include "../include/field.cuh"
#include "../include/params.cuh"
#include <stddef.h>
#include <cuda_runtime.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

// Function declarations
__host__ __device__ int is_zero(const uint64_t *a, size_t field_words);
__host__ __device__ int field_equal(const uint64_t *a, const uint64_t *b, size_t field_words);
__host__ __device__ void field_pow(uint64_t *result, const uint64_t *a, const size_t exponent, size_t field_words);
__host__ __device__ void field_inv(uint64_t *inv, const uint64_t *base, size_t field_words);
__host__ __device__ void field_batch_inverse_and_mul_with_precomputed(uint64_t *result, const uint64_t *elements, const size_t elements_len, const uint64_t *mul, size_t field_words, const uint64_t *precomputed_inverses, size_t num_precomputed);
__host__ __device__ void field_add(uint64_t *sum, const uint64_t *a, const uint64_t *b,size_t field_words);
__host__ __device__ void field_addEqual(uint64_t *a, const uint64_t *b, int field_words);
__host__ __device__ void field_sub(uint64_t *difference, const uint64_t *a, const uint64_t *b,size_t field_words);
__host__ __device__ void field_subEqual(uint64_t *a, const uint64_t *b,size_t field_words);
__host__ __device__ static uint64_t gf64_mul(uint64_t a, uint64_t b);
__host__ __device__ void field_mul(uint64_t *product, const uint64_t *a, const uint64_t *b,size_t field_words);
__host__ __device__ void field_mulEqual(uint64_t *a, const uint64_t *b,size_t field_words);
__host__ __device__ void field_swap_with_tmp(uint64_t *a, uint64_t *b, uint64_t *tmp, const size_t field_bytesize);
__host__ __device__ void field_swap(uint64_t *a, uint64_t *b, const size_t field_bytesize);
__host__ __device__ void field_one(uint64_t *one,size_t field_words);

__constant__ uint64_t zeros[5] = {0, 0, 0, 0, 0};
__constant__ uint64_t one[5] = {0, 0, 0, 0, 1};
__host__ __device__ int device_memcmp(const uint64_t *a, const uint64_t *b, size_t num_elements) {
    for (size_t i = 0; i < num_elements; ++i) {
        if (a[i] != b[i]) {
            return (a[i] < b[i]) ? -1 : 1;  // Mimic the behavior of memcmp
        }
    }
    return 0;
}
// Function definitions
__host__ __device__
int is_zero(const uint64_t *a,size_t field_words)
{
    switch (field_words)
    {
    case 1:
        return device_memcmp(a, zeros, sizeof(uint64_t)) == 0;
    case 2:
        return device_memcmp(a, zeros, sizeof(uint64_t) * 2) == 0;
    case 3:
        return device_memcmp(a, zeros, sizeof(uint64_t) * 3) == 0;
    case 4:
        return device_memcmp(a, zeros, sizeof(uint64_t) * 4) == 0;
    case 5:
        return device_memcmp(a, zeros, sizeof(uint64_t) * 5) == 0;
    default:
        for (size_t i = 0; i < field_words; i++)
        {
            if (a[i] != 0)
                return 0;
        }
        return 1;
    }
}

__host__ __device__ 
int field_equal(const uint64_t *a, const uint64_t *b,size_t field_words)
{
    switch (field_words)
    {
    case 1:
        return !device_memcmp(a, b, sizeof(uint64_t));
    case 2:
        return !device_memcmp(a, b, sizeof(uint64_t) * 2);
    case 3:
        return !device_memcmp(a, b, sizeof(uint64_t) * 3);
    case 4:
        return !device_memcmp(a, b, sizeof(uint64_t) * 4);
    case 5:
        return !device_memcmp(a, b, sizeof(uint64_t) * 5);
    default:
        return !device_memcmp(a, b, sizeof(uint64_t) * field_words);
    }
}

__host__ __device__ 
void field_pow(uint64_t *result, const uint64_t *a, const size_t exponent,size_t field_words)
{
    uint64_t field_one[FIELD_WORDS];
    memset(field_one, 0, (field_words - 1) * sizeof(uint64_t));
    field_one[FIELD_WORDS - 1] = 1;

    memcpy(result, field_one, field_words * sizeof(uint64_t));
    int found_one = 0;

    for (long i = 63; i >= 0; --i)
    {
        if (found_one == 1)
        {
            field_mulEqual(result, result, field_words);
        }

        if (exponent & (1ull << i))
        {
            found_one = 1;
            field_mulEqual(result, a, field_words);
        }
    }
}

__host__ __device__ 
void field_inv(uint64_t *inv, const uint64_t *base,size_t field_words)
{
    if (field_words == 3)
    {
        memset(inv, 0, field_words * sizeof(uint64_t));
        uint64_t *a = (uint64_t *)malloc(field_words * sizeof(uint64_t));
        memcpy(a, base, field_words * sizeof(uint64_t));
        uint64_t *prev_inv = (uint64_t *)malloc(field_words * sizeof(uint64_t));
        for (size_t i = 0; i <= 6; ++i)
        {
            uint64_t *b = (uint64_t *)malloc(field_words * sizeof(uint64_t));
            memcpy(b, a, field_words * sizeof(uint64_t));
            for (size_t j = 0; j < (1ul << i); ++j)
                field_mulEqual(b, b, field_words);
            field_mulEqual(a, b, field_words);

            if (i == 6)
                memcpy(prev_inv, inv, field_words * sizeof(uint64_t));
            if (i == 0)
                memcpy(inv, b, field_words * sizeof(uint64_t));
            else
                field_mulEqual(inv, b, field_words);

            free(b);
        }

        for (size_t i = 0; i < (1ul << 6); ++i)
            field_mulEqual(inv, inv, field_words);
        field_mulEqual(prev_inv, prev_inv, field_words);

        field_mulEqual(inv, prev_inv, field_words);
        field_mulEqual(inv, base, field_words);
        field_mulEqual(inv, base, field_words);

        free(a);
        free(prev_inv);
    }
    else if (field_words == 4)
    {
        // Computing a^{-1} = a^{2^256 - 2}
        memset(inv, 0, field_words * sizeof(uint64_t));
        uint64_t *a = (uint64_t *)malloc(field_words * sizeof(uint64_t));
        memcpy(a, base, field_words * sizeof(uint64_t));
        for (size_t i = 0; i <= 7; ++i)
        {
            /* entering the loop a = base^{2^{2^i}-1} */
            uint64_t *b = (uint64_t *)malloc(field_words * sizeof(uint64_t));
            memcpy(b, a, field_words * sizeof(uint64_t));
            for (size_t j = 0; j < (1ul << i); ++j)
                field_mulEqual(b, b, field_words);
            /* after the loop b = a^{2^{2^i}} = base^{2^{2^i}*(2^{2^i}-1)} */
            field_mulEqual(a, b, field_words);
            /* now a = base^{2^{2^{i+1}}-1} */

            if (i == 0)
                memcpy(inv, b, field_words * sizeof(uint64_t));
            else
                field_mulEqual(inv, b, field_words);
            free(b);
        }

        free(a);
    }
    else if (field_words == 5)
    {
        // Computing a^{-1} = a^{2^320 - 2}
        memset(inv, 0, field_words * sizeof(uint64_t));
        uint64_t *a = (uint64_t *)malloc(field_words * sizeof(uint64_t));
        memcpy(a, base, field_words * sizeof(uint64_t));
        uint64_t *prev_inv = (uint64_t *)malloc(field_words * sizeof(uint64_t));
        for (size_t i = 0; i <= 7; ++i)
        {
            /* entering the loop a = base^{2^{2^i}-1} */
            uint64_t *b = (uint64_t *)malloc(field_words * sizeof(uint64_t));
            memcpy(b, a, field_words * sizeof(uint64_t));
            for (size_t j = 0; j < (1ul << i); ++j)
                field_mulEqual(b, b, field_words);
            /* after the loop b = a^{2^i} = base^{2^{2^i}*(2^{2^i}-1)} */
            field_mulEqual(a, b, field_words);
            /* now a = base^{2^{2^{i+1}}-1} */

            if (i == 6)
                memcpy(prev_inv, inv, field_words * sizeof(uint64_t));
            if (i == 0)
                memcpy(inv, b, field_words * sizeof(uint64_t));
            else
                field_mulEqual(inv, b, field_words);

            free(b);
        }

        /* now inv = base^{2^256-2}, prev_inv = base^{2^64-2} */
        for (size_t i = 0; i < (1ul << 6); ++i)
            field_mulEqual(inv, inv, field_words);
        field_mulEqual(prev_inv, prev_inv, field_words);

        /* now inv = base^{2^320 - 2*2^64}, prev_inv = base^{2*2^64 - 4},
        thus base^{2^320 - 2} = inv * prev_inv * base^{2} */
        field_mulEqual(inv, prev_inv, field_words);
        field_mulEqual(inv, base, field_words);
        field_mulEqual(inv, base, field_words);

        free(a);
        free(prev_inv);
    }
}

__host__ __device__ 
void field_batch_inverse_and_mul_with_precomputed(uint64_t *result, const uint64_t *elements, const size_t elements_len, const uint64_t *mul,size_t field_words, const uint64_t *precomputed_inverses, size_t num_precomputed)
{
    uint64_t c[FIELD_WORDS];
    memcpy(c, &elements[0], field_words * sizeof(uint64_t));
    memcpy(&result[0], c, field_words * sizeof(uint64_t));

    for (size_t i = 1; i < elements_len; ++i)
    {
        field_mulEqual(c, &elements[i * field_words], field_words);
        memcpy(&result[i * field_words], c, field_words * sizeof(uint64_t));
    }

    uint64_t c_inv[FIELD_WORDS];
    if (0 < num_precomputed) {
        memcpy(c_inv, &precomputed_inverses[0 * field_words], field_words * sizeof(uint64_t));
    } else {
        field_inv(c_inv, c, field_words);
    }
    
    field_mulEqual(c_inv, mul, field_words);

    for (size_t i = elements_len - 1; i > 0; --i)
    {
        field_mul(&result[i * field_words], &result[(i - 1) * field_words], c_inv, field_words);
        field_mulEqual(c_inv, &elements[i * field_words], field_words);
    }

    memcpy(&result[0], c_inv, field_words * sizeof(uint64_t));
}

__host__ __device__ 
void field_add(uint64_t *sum, const uint64_t *a, const uint64_t *b,size_t field_words)
{
    for (size_t i = 0; i < field_words; i++)
        sum[i] = a[i] ^ b[i];
}

__host__ __device__ 
void field_addEqual(uint64_t *a, const uint64_t *b, size_t field_words)
{
    for (size_t i = 0; i < field_words; i++)
        a[i] ^= b[i];
}

__host__ __device__ 
void field_sub(uint64_t *difference, const uint64_t *a, const uint64_t *b,size_t field_words)
{
    field_add(difference, a, b, field_words);
}

__host__ __device__ 
void field_subEqual(uint64_t *a, const uint64_t *b,size_t field_words)
{
    field_addEqual(a, b, field_words);
}

__host__ __device__ 
static uint64_t gf64_mul(uint64_t a, uint64_t b)
{
    if (a == 0 || b == 0)
        return 0;
    uint64_t result = a & (-(b & 1));
    for (size_t i = 0; i < 64; i++)
    {
        int ar = ((a >> 63) & 1) * 0x1b;
        a = (a << 1) ^ ar;
        b >>= 1;
        result ^= a & (-(b & 1));
    }
    return result;
}

__host__ __device__ 
void field_mul(uint64_t *product, const uint64_t *a, const uint64_t *b,size_t field_words)
{
    if (field_words == 3)
    {
        product[2] = gf64_mul(a[1], b[0]);
        product[1] = gf64_mul(a[0], b[1]);
        product[2] ^= product[1];
        product[1] = product[2];
        product[0] = gf64_mul(a[0], b[0]);
        product[1] ^= product[0];

        product[2] ^= gf64_mul(a[2], b[2]);
        product[1] ^= gf64_mul(a[2], b[1]);
        product[1] ^= gf64_mul(a[1], b[2]);
        product[0] ^= gf64_mul(a[2], b[0]);
        product[0] ^= gf64_mul(a[1], b[1]);
        product[0] ^= gf64_mul(a[0], b[2]);
    }
    else if (field_words == 4)
    {
        product[0] = gf64_mul(a[0], b[0]);
        product[1] = gf64_mul(a[0], b[1]);
        product[1] ^= gf64_mul(a[1], b[0]);
        product[2] = gf64_mul(a[0], b[2]);
        product[2] ^= gf64_mul(a[1], b[1]);
        product[2] ^= gf64_mul(a[2], b[0]);
        product[3] = product[2];
        product[2] ^= product[1];
        product[1] ^= product[0];

        product[3] ^= gf64_mul(a[3], b[3]);
        product[2] ^= gf64_mul(a[3], b[2]);
        product[2] ^= gf64_mul(a[2], b[3]);
        product[1] ^= gf64_mul(a[3], b[1]);
        product[1] ^= gf64_mul(a[2], b[2]);
        product[1] ^= gf64_mul(a[1], b[3]);
        product[0] ^= gf64_mul(a[3], b[0]);
        product[0] ^= gf64_mul(a[2], b[1]);
        product[0] ^= gf64_mul(a[1], b[2]);
        product[0] ^= gf64_mul(a[0], b[3]);
    }
    else if (field_words == 5)
    {
        product[4] = gf64_mul(a[0], b[3]);
        product[4] ^= gf64_mul(a[1], b[2]);
        product[4] ^= gf64_mul(a[2], b[1]);
        product[4] ^= gf64_mul(a[3], b[0]);
        product[3] = gf64_mul(a[0], b[2]);
        product[3] ^= gf64_mul(a[1], b[1]);
        product[3] ^= gf64_mul(a[2], b[0]);
        product[2] = gf64_mul(a[0], b[1]);
        product[2] ^= gf64_mul(a[1], b[0]);
        product[1] = gf64_mul(a[0], b[0]);
        product[0] = product[2];
        product[2] ^= product[4];
        product[4] ^= product[1];
        product[2] ^= product[1];
        product[1] ^= product[3];

        product[4] ^= gf64_mul(a[4], b[4]);
        product[3] ^= gf64_mul(a[3], b[4]);
        product[3] ^= gf64_mul(a[4], b[3]);
        product[2] ^= gf64_mul(a[2], b[4]);
        product[2] ^= gf64_mul(a[3], b[3]);
        product[2] ^= gf64_mul(a[4], b[2]);
        product[1] ^= gf64_mul(a[1], b[4]);
        product[1] ^= gf64_mul(a[2], b[3]);
        product[1] ^= gf64_mul(a[3], b[2]);
        product[1] ^= gf64_mul(a[4], b[1]);
        product[0] ^= gf64_mul(a[0], b[4]);
        product[0] ^= gf64_mul(a[1], b[3]);
        product[0] ^= gf64_mul(a[2], b[2]);
        product[0] ^= gf64_mul(a[3], b[1]);
        product[0] ^= gf64_mul(a[4], b[0]);
    }
}

__host__ __device__ 
void field_mulEqual(uint64_t *a, const uint64_t *b,size_t field_words)
{
    uint64_t *temp;
    uint64_t temp_stack[5];
    if (field_words <= 5)
    {
        temp = (uint64_t *)&temp_stack;
    }
    else
    {
        temp = (uint64_t *)malloc(field_words * sizeof(uint64_t));
    }
    field_mul(temp, a, b, field_words);
    memcpy(a, temp, field_words * sizeof(uint64_t));
    if (field_words > 5)
    {
        free(temp);
    }
}

__host__ __device__ 
void field_swap_with_tmp(uint64_t *a, uint64_t *b, uint64_t *tmp, const size_t field_bytesize)
{
    memcpy(tmp, a, field_bytesize);
    memcpy(a, b, field_bytesize);
    memcpy(b, tmp, field_bytesize);
}

__host__ __device__ 
void field_swap(uint64_t *a, uint64_t *b, const size_t field_bytesize)
{
    uint64_t *tmp = (uint64_t *)malloc(field_bytesize);
    field_swap_with_tmp(a, b, tmp, field_bytesize);
    free(tmp);
}

__host__ __device__ 
void field_one(uint64_t *one,size_t field_words)
{
    memset(one, 0, (field_words - 1) * sizeof(uint64_t));
    one[field_words - 1] = 1;
}

__host__ __device__ 
void field_neg(uint64_t *result, const uint64_t *a,size_t field_words)
{
    memset(result, 0, field_words * sizeof(uint64_t));
    field_sub(result, result, a, field_words);
}

__host__ __device__ 
void all_subset_sums(uint64_t *sums, const uint64_t *elements, int num_elements, const uint64_t *shift,size_t field_words)
{
    size_t num_sums = 1 << num_elements;
    memset(sums, 0, num_sums * field_words * sizeof(uint64_t));
    for (size_t i = 0; i < num_sums; i++)
    {
        for (size_t j = 0; j < num_elements; j++)
        {
            if (i & (1 << j))
            {
                field_addEqual(&sums[i * field_words], &elements[j * field_words], field_words);
            }
        }
        field_addEqual(&sums[i * field_words], shift, field_words);
    }
}

   
   // Test function - keep commented but do not delete
   // int main() {
   //     int /*replaced from size_t*/ field_words = 3;
   
   //     uint64_t a[3] = {0x0000000000002001, 0x0, 0x0};
   //     uint64_t b[3] = {0x4000000000000000, 0x0, 0x1};
   
   //     printf("Field element a: ");
   //     for (int /*replaced from size_t*/ i = 0; i < field_words; i++) {
   //         printf("%016llx ", a[i]);
   //     }
   //     printf("\n");
   
   //     printf("Field element b: ");
   //     for (int /*replaced from size_t*/ i = 0; i < field_words; i++) {
   //         printf("%016llx ", b[i]);
   //     }
   //     printf("\n");
   
   //     uint64_t result1[3], result2[3];
   //     field_mul(result1, a, a, field_words);
   //     field_mul(result1, result1, result1, field_words);
   //     field_pow(result2, a, 4, field_words);
   
   //     printf("a*a*a*a == a**4: %s\n", field_equal(result1, result2, field_words) ? "true" : "false");
   
   //     uint64_t inv_a[3];
   //     field_inv(inv_a, a, field_words);
   //     printf("Inverse of a (inv_a): ");
   //     for (int /*replaced from size_t*/ i = 0; i < field_words; i++) {
   //         printf("%016llx ", inv_a[i]);
   //     }
   //     printf("\n");
   
   //     // a.inverse() * a == 1
   //     field_mul(result1, inv_a, a, field_words);
   //     printf("a.inverse() * a: ");
   //     for (int /*replaced from size_t*/ i = 0; i < field_words; i++) {
   //         printf("%016llx ", result1[i]);
   //     }
   
   //     uint64_t inv_b[3], div_ab[3], mul_abinv[3];
   //     field_inv(inv_b, b, field_words);
   //     field_mul(div_ab, a, inv_b, field_words);
   
   //     field_mulEqual(result1, b, field_words);
   //     field_inv(inv_b, b, field_words);
   //     field_mul(mul_abinv, a, inv_b, field_words);
   
   //     printf("a / b == a * b.inverse(): %s\n", field_equal(div_ab, mul_abinv, field_words) ? "true" : "false");
   
   //     printf("a + b: ");
   //     field_add(result1, a, b, field_words);
   //     for (int /*replaced from size_t*/ i = 0; i < field_words; i++) {
   //         printf("%016llx ", result1[i]);
   //     }
   //     printf("\n");
   
   //     printf("a - b: ");
   //     field_sub(result1, a, b, field_words);
   //     for (int /*replaced from size_t*/ i = 0; i < field_words; i++) {
   //         printf("%016llx ", result1[i]);
   //     }
   //     printf("\n");
   
   //     printf("a * b: ");
   //     field_mul(result1, a, b, field_words);
   //     for (int /*replaced from size_t*/ i = 0; i < field_words; i++) {
   //         printf("%016llx ", result1[i]);
   //     }
   //     printf("\n");
   
   //     printf("a / b: ");
   //     for (int /*replaced from size_t*/ i = 0; i < field_words; i++) {
   //         printf("%016llx ", div_ab[i]);
   //     }
   //     printf("\n");
   
   //     field_pow(result1, a, 5, field_words);
   //     printf("a**5: ");
   //     for (int /*replaced from size_t*/ i = 0; i < field_words; i++) {
   //         printf("%016llx ", result1[i]);
   //     }
   //     printf("\n");
   
   //     uint64_t zero[3] = {0, 0, 0};
   //     printf("Zero check (field zero): ");
   //     for (int /*replaced from size_t*/ i = 0; i < field_words; i++) {
   //         printf("%016llx ", zero[i]);
   //     }
   //     printf("\n");
   
   //     return 0;
   // }