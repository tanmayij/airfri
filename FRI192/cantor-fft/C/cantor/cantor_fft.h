#ifndef _CANTOR_FFT_H_
#define _CANTOR_FFT_H_
#include <immintrin.h>
    __m128i* cantor_fft_gf2128(__m128i* poly, unsigned n_term);

    __m128i* cantor_fft_hc_circuits_gf2128(__m128i* poly, unsigned n_term); //hardcodes the circuits
    __m128i* cantor_fft_hc_cache_gf2128(__m128i* poly, unsigned n_term); //hardcodes the circuits
    
#endif
