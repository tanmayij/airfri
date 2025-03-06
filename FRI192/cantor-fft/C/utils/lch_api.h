#ifndef _LCH_API_H_
#define _LCH_API_H_
#include <immintrin.h>
    /*
    The bitpolymul doesn't have butterfly_net_clmul function, since it is for multiplying polynomials. 
    Therefore, they have butterfly_net_half_inp_clmul instead.
    */
    __m128i* lch_fft_gf2128(__m128i* fx, unsigned n_term);
#endif
