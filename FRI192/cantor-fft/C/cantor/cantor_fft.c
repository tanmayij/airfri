#include "cantor_fft.h"
#include <stdio.h>
#include <string.h>    // Header for memcpy

#include "s_i.h"
#include "utils/utils.h"
#include "bitpolymul/bitmat_prod.h"
#include "bitpolymul/gf2128_cantor_iso.h"
#include "bitpolymul/gfext_aesni.h"



#define LOG2(X) ((unsigned) (8*sizeof (unsigned long long) - __builtin_clzll((X)) - 1))

__m128i* cantor_fft_gf2128(__m128i* fx, unsigned n_term){
    #ifdef NCOPY_POLY //For small polynomials copy to a new array to keep the input values unchanged
    __m128i* poly = fx;
    #else
    __m128i* poly = (__m128i*)aligned_alloc( 32 , sizeof(__m128i)*n_term );
    poly = memcpy( poly, fx, sizeof(__m128i)*n_term );
    #endif
    unsigned m = LOG2(n_term);
    unsigned S_index = m, input_size = n_term, n_modules = 1;
    unsigned j, t, offset, offset2, half_input_size, half_half_input_size;
    __m128i mult_factor, poly_k;
    __m256i* poly256;
    // unsigned nz_S[15];
    const unsigned* nz_S;
    for (unsigned r = 0; r < m-1; ++r){
        // Computing S_i
        S_index--; 
        nz_S = s_i[S_index]; t=n_terms[S_index];
        offset = 0;
        half_input_size = input_size >> 1;
        half_half_input_size = half_input_size>>1;
        for (unsigned module = 0; module < n_modules; ++module){
            mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
            offset2 = offset + half_input_size;
            for (unsigned k = offset2 + half_input_size - 1; k >= offset2; --k){
                poly_k = poly[k];
                for (unsigned i = 0; i < t; ++i){
                    poly[k - nz_S[i]] ^= poly_k;
                }
                poly[k-half_input_size] ^= _gf2ext128_mul_sse(poly_k, mult_factor);
            }
            poly256 =(__m256i*) (poly+offset);
            for (j = 0; j < half_half_input_size; ++j) { // we use half_half_input_size since the steps are 256-bits 
                // poly256[j + half_half_input_size] = poly256[j] ^= _gf2ext128_mul_2x1_avx2( poly256[j + half_half_input_size] , mult_factor );
                poly256[j + half_half_input_size] ^= poly256[j];  
            }      
            offset += input_size;
        }
        input_size = half_input_size;
        n_modules <<= 1;
    }
    
    offset = 0;
    // r = m (last layer)
    for (unsigned module = 0; module < n_modules; ++module){
        mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
        offset2 = offset + half_input_size;
        poly[offset] ^= _gf2ext128_mul_sse(poly[offset+1], mult_factor);
        poly[offset+1] ^= poly[offset];            
        offset += 2;
    }

    return poly;
}

#include "div_by_s_i.h"
__m128i* cantor_fft_hc_circuits_gf2128(__m128i* fx, unsigned n_term){
    #ifdef NCOPY_POLY //For small polynomials copy to a new array to keep the input values unchanged
    __m128i* poly = fx;
    #else
    __m128i* poly = (__m128i*)aligned_alloc( 32 , sizeof(__m128i)*n_term );
    poly = memcpy( poly, fx, sizeof(__m128i)*n_term );
    #endif
    unsigned m = LOG2(n_term);
    unsigned S_index = m, input_size = n_term, n_modules = 1;
    unsigned j, t, offset, offset2, half_input_size, half_half_input_size;
    __m128i mult_factor, poly_k;
    __m128i* poly_offset;
    __m256i* poly256;
    // unsigned nz_S[15];
    const unsigned* nz_S;
    for (unsigned r = 0; r < m-5; ++r){
        // Computing S_i
        S_index--; 
        nz_S = s_i[S_index]; t=n_terms[S_index];
        offset = 0;
        half_input_size = input_size >> 1;
        half_half_input_size = half_input_size>>1;
        for (unsigned module = 0; module < n_modules; ++module){
            mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
            offset2 = offset + half_input_size;
            for (unsigned k = offset2 + half_input_size - 1; k >= offset2; --k){
                poly_k = poly[k];
                for (unsigned i = 0; i < t; ++i){
                    poly[k - nz_S[i]] ^= poly_k;
                }
                poly[k-half_input_size] ^= _gf2ext128_mul_sse(poly_k, mult_factor);
            }
            poly256 =(__m256i*) (poly+offset);
            for (j = 0; j < half_half_input_size; ++j) { // we use half_half_input_size since the steps are 256-bits 
                // poly256[j + half_half_input_size] = poly256[j] ^= _gf2ext128_mul_2x1_avx2( poly256[j + half_half_input_size] , mult_factor );
                poly256[j + half_half_input_size] ^= poly256[j];  
            }      
            offset += input_size;
        }
        input_size = half_input_size;
        n_modules <<= 1;
    }

    // r = m-4 (S_4) (input_size = 32)
    offset = 0;
    half_input_size = input_size >> 1;
    for (unsigned module = 0; module < n_modules; ++module){    
        mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
        poly_offset = (__m128i*) (poly + offset);  
        DIV_S4(mult_factor);
        offset += 32;
    }
    n_modules <<= 1;

    // r = m-3 (S_3) (input_size = 16)
    offset = 0;
    half_input_size = input_size >> 1;
    for (unsigned module = 0; module < n_modules; ++module){    
        mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
        poly_offset = (__m128i*) (poly + offset);  
        DIV_S3(mult_factor);
        offset += 16;
    }
    n_modules <<= 1;

    // r = m-2 (S_2) (input_size = 8)
    offset = 0;
    half_input_size = input_size >> 1;
    for (unsigned module = 0; module < n_modules; ++module){    
        mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
        poly_offset = (__m128i*) (poly + offset);  
        DIV_S2(mult_factor);
        offset += 8;
    }
    n_modules <<= 1;

    // r = m-1 (S_1) (input_size = 4)
    offset = 0;
    half_input_size = input_size >> 1;
    for (unsigned module = 0; module < n_modules; ++module){    
        mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
        poly_offset = (__m128i*) (poly + offset);                       
        DIV_S1(mult_factor);
        offset += 4;
    }
    n_modules <<= 1;

    // r = m (last layer) (S_0)
    offset = 0;
    for (unsigned module = 0; module < n_modules; ++module){    
        mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
        poly_offset = (__m128i*) (poly + offset);                       
        DIV_S0(mult_factor); 
        offset += 2;
    }

    return poly;
}



__m128i* cantor_fft_hc_cache_gf2128(__m128i* fx, unsigned n_term){
    #ifdef NCOPY_POLY //For small polynomials copy to a new array to keep the input values unchanged
    __m128i* poly = fx;
    #else
    __m128i* poly = (__m128i*)aligned_alloc( 32 , sizeof(__m128i)*n_term );
    poly = memcpy( poly, fx, sizeof(__m128i)*n_term );
    #endif
    unsigned m = LOG2(n_term);
    unsigned S_index = m, input_size = n_term, n_modules = 1;
    unsigned j, t, offset, offset2, half_input_size, half_half_input_size;
    __m128i mult_factor, poly_k;
    __m128i* poly_offset;
    __m256i* poly256;
    // unsigned nz_S[15];
    const unsigned* nz_S;
    for (unsigned r = 0; r < m-5; ++r){
        // Computing S_i
        S_index--; 
        nz_S = s_i[S_index]; t=n_terms[S_index];
        offset = 0;
        half_input_size = input_size >> 1;
        half_half_input_size = half_input_size>>1;
        for (unsigned module = 0; module < n_modules; ++module){
            mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
            offset2 = offset + half_input_size;
            for (unsigned k = offset2 + half_input_size - 1; k >= offset2; --k){
                poly_k = poly[k];
                for (unsigned i = 0; i < t; ++i){
                    poly[k - nz_S[i]] ^= poly_k;
                }
                poly[k-half_input_size] ^= _gf2ext128_mul_sse(poly_k, mult_factor);
            }
            poly256 =(__m256i*) (poly+offset);
            for (j = 0; j < half_half_input_size; ++j) { // we use half_half_input_size since the steps are 256-bits 
                // poly256[j + half_half_input_size] = poly256[j] ^= _gf2ext128_mul_2x1_avx2( poly256[j + half_half_input_size] , mult_factor );
                poly256[j + half_half_input_size] ^= poly256[j];  
            }      
            offset += input_size;
        }
        input_size = half_input_size;
        n_modules <<= 1;
    }

    // r = m-4 (S_4) (input_size = 32)
    offset = 0;
    half_input_size = input_size >> 1;
    for (unsigned module = 0; module < n_modules; ++module){    
        mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
        poly_offset = (__m128i*) (poly + offset);  
        DIV_S4(mult_factor);
        offset += 32;
    }
    n_modules <<= 1;

    // r = m-3 (S_3) (input_size = 16)
    offset = 0;
    half_input_size = input_size >> 1;
    for (unsigned module = 0; module < n_modules; ++module){    
        mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
        poly_offset = (__m128i*) (poly + offset);  
        DIV_S3(mult_factor);
        offset += 16;
    }
    n_modules <<= 1;

    // r = m-2 (S_2) (input_size = 8)
    offset = 0;
    half_input_size = input_size >> 1;
    for (unsigned module = 0; module < n_modules; ++module){    
        mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
        poly_offset = (__m128i*) (poly + offset);  
        DIV_S2(mult_factor);
        offset += 8;
    }
    n_modules <<= 1;

    // r = m-1 (S_1) (input_size = 4)
    offset = 0;
    half_input_size = input_size >> 1;
    for (unsigned module = 0; module < n_modules; ++module){    
        // unsigned module_2 = module<<1;
        // unsigned module_4 = module<<2;
        mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
        poly_offset = (__m128i*) (poly + offset);                       
        DIV_S1(mult_factor);
        
        mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<2);
        DIV_S0(mult_factor); 
        poly_offset += 2;                       
        // xmm_dump(&mult_factor, 1);
        // mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, (module<<2)+2);
        mult_factor ^= _mm_load_si128((const __m128i*) gfCantorto2128_8R + 2);
        // xmm_dump(&mult_factor, 1);
        // printf("\n");
        // printf("\n");

        // mult_factor ^= _mm_set_epi32(0, 0, 0, 1);
        DIV_S0(mult_factor); 
        offset += 4;
    }
    // n_modules <<= 1;

    // r = m (last layer) (S_0)
    // offset = 0;
    // for (unsigned module = 0; module < n_modules; ++module){    
    //     mult_factor = bitmat_prod_accu_64x128_M8R_sse(_mm_setzero_si128(),  gfCantorto2128_8R, module<<1);
    //     poly_offset = (__m128i*) (poly + offset);                       
    //     DIV_S0(mult_factor); 
    //     offset += 2;
    // }

    return poly;
}