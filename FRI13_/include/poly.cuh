#ifndef POLY_H
#define POLY_H

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include <stdbool.h>

#include "domain.cuh"
#include "field.cuh"
//#include "fft.h"
//#include "fri.cuh"

// Polynomial operations
extern void fft(uint64_t *result, const uint64_t *v, const size_t v_len, const size_t field_words, const size_t domain_size, const uint64_t *domain_shift);
int poly_add(uint64_t *sum, const uint64_t *a, const size_t a_len, const uint64_t *b, const size_t b_len, const size_t field_words);

void poly_addEqual(uint64_t *a, const size_t a_len, const uint64_t *b, const size_t b_len, const size_t field_words);

void poly_scalarMul(uint64_t *product, const uint64_t *a, const size_t a_len, const uint64_t *b, const size_t field_words);

void poly_scalarMulEqual(uint64_t *a, const size_t a_len, const uint64_t *b, const size_t field_words);

void poly_mul(uint64_t *product, const uint64_t *a, const size_t a_len, const uint64_t *b, const size_t b_len, const size_t field_words);

void poly_mulEqual(uint64_t *a, const size_t a_len, const uint64_t *b, const size_t b_len, const size_t field_words);

void poly_div(uint64_t *quotient, uint64_t *remainder, const uint64_t *a, const size_t a_len, const uint64_t *b, const size_t b_len, const size_t field_words);

void poly_eval(uint64_t *result, const uint64_t *poly, const size_t poly_len, const uint64_t *x, const size_t field_words);

void poly_eval_over_domain(uint64_t *result, const uint64_t *poly, const size_t poly_len, const Domain *domain, const size_t field_words, const size_t field_bytesize);

int poly_deg(const uint64_t *poly, size_t poly_len, size_t field_words);

bool is_zero_poly(const uint64_t *poly, const size_t poly_len, const size_t field_words);

//Polynomial* interpolate_domain_polynomial(uint64_t **domain, uint64_t **values, int num_points, size_t field_words);
void interpolate_domain_single(uint64_t *result, uint64_t **domain, uint64_t *values, size_t domain_length, size_t field_words);
//bool test_colinearity(uint64_t **points, size_t num_points, size_t field_words);

void shift_poly(uint64_t *result, const uint64_t *poly, size_t poly_len, const uint64_t *offset, size_t field_words);

void scale_poly(uint64_t *result, const uint64_t *poly, size_t poly_len, const uint64_t *factor, size_t field_words);

#endif // POLY_H