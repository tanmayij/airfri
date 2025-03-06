#ifndef _UTILS_H_
#define _UTILS_H_
#include <immintrin.h>
#include <smmintrin.h>
#include "bitpolymul/byte_inline_func.h"

static inline bool are_equal_128bit(__m128i in1, __m128i in2) {in1 ^= in2; return _mm_testz_si128(in1, in1);}
static inline bool are_equal_vec_128bit(__m128i* in1, __m128i* in2, unsigned n) {bool eq=true;for(unsigned i=0;i<n;++i) eq&=are_equal_128bit(in1[i],in2[i]);return eq;}
static inline void print_vector_128bit_name(__m128i * vec, unsigned n, const char * vec_name){printf("%s :", vec_name); xmm_dump(vec, n); puts("");}
static inline void print_vector_128bit(__m128i * vec, unsigned n){print_vector_128bit_name(vec, n, "vec");}

__m128i* naive_evaluate(__m128i* fx, unsigned n_term);

__m128i* random_polynomial_gf2128(unsigned n);
__m128i* zero_vec_128bit(unsigned n);



bool validate_cantor_basis();

#endif
