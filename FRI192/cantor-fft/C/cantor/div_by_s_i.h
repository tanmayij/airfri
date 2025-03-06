#ifndef _DIV_BY_S_I_H_
#define _DIV_BY_S_I_H_

#define DIV_S0(mult_factor) do {                                \
    poly_offset[0] ^= _gf2ext128_mul_sse(poly_offset[1], (mult_factor));  \
    poly_offset[1] ^= poly_offset[0];                                    \
} while (0)

#define DIV_S1(mult_factor) do {                                \
        poly_offset[2] ^= poly_offset[3];                               \
        poly_offset[1] ^= _gf2ext128_mul_sse(poly_offset[3], mult_factor) ^ poly_offset[2];                                         \
        poly_offset[0] ^= _gf2ext128_mul_sse(poly_offset[2], mult_factor);                                                          \
        poly_offset[2] ^= poly_offset[0]; \
        poly_offset[3] ^= poly_offset[1]; \
} while (0)

#define DIV_S2(mult_factor) do {                                \
        poly_offset[4] ^= poly_offset[7];                       \
        poly_offset[3] ^= poly_offset[6];                       \
        poly_offset[2] ^= poly_offset[5];                       \
        poly_offset[1] ^= poly_offset[4];                       \
        poly_offset[3] ^= _gf2ext128_mul_sse(poly_offset[7], mult_factor);  \
        poly_offset[2] ^= _gf2ext128_mul_sse(poly_offset[6], mult_factor);  \
        poly_offset[1] ^= _gf2ext128_mul_sse(poly_offset[5], mult_factor);  \
        poly_offset[0] ^= _gf2ext128_mul_sse(poly_offset[4], mult_factor);  \
        poly_offset[4] ^= poly_offset[0];                       \
        poly_offset[5] ^= poly_offset[1];                       \
        poly_offset[6] ^= poly_offset[2];                       \
        poly_offset[7] ^= poly_offset[3];                       \
} while (0)


#define DIV_S3(mult_factor) do {                                \
        poly_offset[8]  ^= poly_offset[15];                       \
        poly_offset[9]  ^= poly_offset[15];                       \
        poly_offset[11] ^= poly_offset[15];                       \
        poly_offset[7]  ^= poly_offset[14];                       \
        poly_offset[8]  ^= poly_offset[14];                       \
        poly_offset[10] ^= poly_offset[14];                       \
        poly_offset[6]  ^= poly_offset[13];                       \
        poly_offset[7]  ^= poly_offset[13];                       \
        poly_offset[9]  ^= poly_offset[13];                       \
        poly_offset[5]  ^= poly_offset[12];                       \
        poly_offset[6]  ^= poly_offset[12];                       \
        poly_offset[8]  ^= poly_offset[12];                       \
        poly_offset[4]  ^= poly_offset[11];                       \
        poly_offset[5]  ^= poly_offset[11];                       \
        poly_offset[7]  ^= poly_offset[11];                       \
        poly_offset[3]  ^= poly_offset[10];                       \
        poly_offset[4]  ^= poly_offset[10];                       \
        poly_offset[6]  ^= poly_offset[10];                       \
        poly_offset[2]  ^= poly_offset[9];                       \
        poly_offset[3]  ^= poly_offset[9];                       \
        poly_offset[5]  ^= poly_offset[9];                       \
        poly_offset[1]  ^= poly_offset[8];                       \
        poly_offset[2]  ^= poly_offset[8];                       \
        poly_offset[4]  ^= poly_offset[8];                       \
        poly_offset[7] ^= _gf2ext128_mul_sse(poly_offset[15], mult_factor);  \
        poly_offset[6] ^= _gf2ext128_mul_sse(poly_offset[14], mult_factor);  \
        poly_offset[5] ^= _gf2ext128_mul_sse(poly_offset[13], mult_factor);  \
        poly_offset[4] ^= _gf2ext128_mul_sse(poly_offset[12], mult_factor);  \
        poly_offset[3] ^= _gf2ext128_mul_sse(poly_offset[11], mult_factor);  \
        poly_offset[2] ^= _gf2ext128_mul_sse(poly_offset[10], mult_factor);  \
        poly_offset[1] ^= _gf2ext128_mul_sse(poly_offset[9], mult_factor);   \
        poly_offset[0] ^= _gf2ext128_mul_sse(poly_offset[8], mult_factor);   \
        poly_offset[8]  ^= poly_offset[0];                       \
        poly_offset[9]  ^= poly_offset[1];                       \
        poly_offset[10] ^= poly_offset[2];                       \
        poly_offset[11] ^= poly_offset[3];                       \
        poly_offset[12] ^= poly_offset[4];                       \
        poly_offset[13] ^= poly_offset[5];                       \
        poly_offset[14] ^= poly_offset[6];                       \
        poly_offset[15] ^= poly_offset[7];                       \
} while (0)


#define DIV_S4(mult_factor) do {                                \
        poly_offset[16] ^= poly_offset[31];                       \
        poly_offset[15] ^= poly_offset[30];                       \
        poly_offset[14] ^= poly_offset[29];                       \
        poly_offset[13] ^= poly_offset[28];                       \
        poly_offset[12] ^= poly_offset[27];                       \
        poly_offset[11] ^= poly_offset[26];                       \
        poly_offset[10] ^= poly_offset[25];                       \
        poly_offset[9]  ^= poly_offset[24];                       \
        poly_offset[8]  ^= poly_offset[23];                       \
        poly_offset[7]  ^= poly_offset[22];                       \
        poly_offset[6]  ^= poly_offset[21];                       \
        poly_offset[5]  ^= poly_offset[20];                       \
        poly_offset[4]  ^= poly_offset[19];                       \
        poly_offset[3]  ^= poly_offset[18];                       \
        poly_offset[2]  ^= poly_offset[17];                       \
        poly_offset[1]  ^= poly_offset[16];                       \
        poly_offset[15] ^= _gf2ext128_mul_sse(poly_offset[31], mult_factor);  \
        poly_offset[14] ^= _gf2ext128_mul_sse(poly_offset[30], mult_factor);  \
        poly_offset[13] ^= _gf2ext128_mul_sse(poly_offset[29], mult_factor);  \
        poly_offset[12] ^= _gf2ext128_mul_sse(poly_offset[28], mult_factor);  \
        poly_offset[11] ^= _gf2ext128_mul_sse(poly_offset[27], mult_factor);  \
        poly_offset[10] ^= _gf2ext128_mul_sse(poly_offset[26], mult_factor);  \
        poly_offset[9]  ^= _gf2ext128_mul_sse(poly_offset[25], mult_factor);  \
        poly_offset[8]  ^= _gf2ext128_mul_sse(poly_offset[24], mult_factor);  \
        poly_offset[7]  ^= _gf2ext128_mul_sse(poly_offset[23], mult_factor);  \
        poly_offset[6]  ^= _gf2ext128_mul_sse(poly_offset[22], mult_factor);  \
        poly_offset[5]  ^= _gf2ext128_mul_sse(poly_offset[21], mult_factor);  \
        poly_offset[4]  ^= _gf2ext128_mul_sse(poly_offset[20], mult_factor);  \
        poly_offset[3]  ^= _gf2ext128_mul_sse(poly_offset[19], mult_factor);  \
        poly_offset[2]  ^= _gf2ext128_mul_sse(poly_offset[18], mult_factor);  \
        poly_offset[1]  ^= _gf2ext128_mul_sse(poly_offset[17], mult_factor);  \
        poly_offset[0]  ^= _gf2ext128_mul_sse(poly_offset[16], mult_factor);  \
        poly_offset[16] ^= poly_offset[0];                       \
        poly_offset[17] ^= poly_offset[1];                       \
        poly_offset[18] ^= poly_offset[2];                       \
        poly_offset[19] ^= poly_offset[3];                       \
        poly_offset[20] ^= poly_offset[4];                       \
        poly_offset[21] ^= poly_offset[5];                       \
        poly_offset[22] ^= poly_offset[6];                       \
        poly_offset[23] ^= poly_offset[7];                       \
        poly_offset[24] ^= poly_offset[8];                       \
        poly_offset[25] ^= poly_offset[9];                       \
        poly_offset[26] ^= poly_offset[10];                      \
        poly_offset[27] ^= poly_offset[11];                      \
        poly_offset[28] ^= poly_offset[12];                      \
        poly_offset[29] ^= poly_offset[13];                      \
        poly_offset[30] ^= poly_offset[14];                      \
        poly_offset[31] ^= poly_offset[15];                      \
} while (0)



#endif
