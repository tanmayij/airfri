#include "utils.h"


#include "bitpolymul/gf2128_cantor_iso.h"
#include "bitpolymul/gfext_aesni.h"
#include "bitpolymul/bitmat_prod.h"
#include <stdio.h>

__m128i* naive_evaluate(__m128i* fx, unsigned n_term){
    __m128i* evaluations = zero_vec_128bit(n_term);
    for (unsigned i=0; i<n_term; i++){
         __m128i p = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, i);
        for (unsigned j=n_term; j--;){
            evaluations[i] = _gf2ext128_mul_sse(evaluations[i], p);
            evaluations[i] ^= fx[j];
        }
    }
    return evaluations;
}

__m128i* random_polynomial_gf2128(unsigned n){
    __m128i* fx = (__m128i*)aligned_alloc( 32 , sizeof(__m128i)*n );
	if( NULL == fx ) { printf("alloc fail.\n"); exit(-1); }
	for (unsigned i = 0; i < n; ++i) fx[i] = _mm_set_epi32(rand(), rand(), rand(), rand());
    return fx;
} 

__m128i* zero_vec_128bit(unsigned n){
    __m128i* vec = (__m128i*)aligned_alloc( 32 , sizeof(__m128i)*n );
	if( NULL == vec ) { printf("alloc fail.\n"); exit(-1); }
	for (unsigned i = 0; i < n; ++i) vec[i] = _mm_setzero_si128();
    return vec;
} 


bool validate_cantor_basis(){
    bool is_cantor_basis = true;
    __m128i prev_beta = _mm_set_epi32(0, 0, 0, 1);
    for (int k = 0; k < 4; ++k){
        for (int i = 0; i < 8; ++i){
            int j = 1<<i;
            __m128i* beta = (__m128i *) gfCantorto2128_8R+j + k*256;
            __m128i expected_prev_beta = _gf2ext128_mul_sse(beta[0], beta[0]);
            expected_prev_beta ^= beta[0];
            is_cantor_basis &= are_equal_128bit(expected_prev_beta, prev_beta);    
            prev_beta = *beta;
        }
    }
    return is_cantor_basis;
}