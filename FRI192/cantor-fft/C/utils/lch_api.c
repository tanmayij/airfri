#include "lch_api.h"
#include "utils.h"
#include <stdio.h>
#include "bitpolymul/bc.h"
#include "bitpolymul/butterfly_net.h"

#include <string.h>    // Header for memcpy


__m128i* lch_fft_gf2128(__m128i* fx, unsigned n_term){
    #ifdef NCOPY_POLY //For small polynomials copy to a new array to keep the input values unchanged
    __m128i* poly = fx;
    #else
    __m128i* poly = (__m128i*)aligned_alloc( 32 , sizeof(__m128i)*n_term );
    poly = memcpy( poly, fx, sizeof(__m128i)*n_term );
    #endif
    bc_to_lch_128( (bc_sto_t*) poly,  n_term );
    butterfly_net_clmul( poly, n_term );

    return poly;
}