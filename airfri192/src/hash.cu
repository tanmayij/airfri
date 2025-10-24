/* -------------------------------------------------------------------------
 * Works when compiled for either 32-bit or 64-bit targets, optimized for 
 * 64 bit.
 *
 * Canonical implementation of Init/Update/Finalize for SHA-3 byte input. 
 *
 * SHA3-256, SHA3-384, SHA-512 are implemented. SHA-224 can easily be added.
 *
 * Based on code from http://keccak.noekeon.org/ .
 *
 * I place the code that I wrote into public domain, free to use. 
 *
 * I would appreciate if you give credits to this work if you used it to 
 * write or test * your code.
 *
 * Aug 2015. Andrey Jivsov. crypto@brainhub.org
 * ---------------------------------------------------------------------- */
 #include <cstdio>
 #include <cstdint>
 #include <cstring>
 
 #include "../include/hash.cuh"
 
 #define SHA3_ASSERT( x )
 #define SHA3_TRACE( format, ...)
 #define SHA3_TRACE_BUF(format, buf, l)
 
 /* 
  * This flag is used to configure "pure" Keccak, as opposed to NIST SHA3.
  */
 #define SHA3_USE_KECCAK_FLAG 0x80000000
 #define SHA3_CW(x) ((x) & (~SHA3_USE_KECCAK_FLAG))
 
 #if defined(_MSC_VER)
 #define SHA3_CONST(x) x
 #else
 #define SHA3_CONST(x) x##L
 #endif
 
 #ifndef SHA3_ROTL64
 #define SHA3_ROTL64(x, y) \
     (((x) << (y)) | ((x) >> ((sizeof(uint64_t)*8) - (y))))
 #endif
 
 __device__ static const uint64_t keccakf_rndc[24] = {
     SHA3_CONST(0x0000000000000001UL), SHA3_CONST(0x0000000000008082UL),
     SHA3_CONST(0x800000000000808aUL), SHA3_CONST(0x8000000080008000UL),
     SHA3_CONST(0x000000000000808bUL), SHA3_CONST(0x0000000080000001UL),
     SHA3_CONST(0x8000000080008081UL), SHA3_CONST(0x8000000000008009UL),
     SHA3_CONST(0x000000000000008aUL), SHA3_CONST(0x0000000000000088UL),
     SHA3_CONST(0x0000000080008009UL), SHA3_CONST(0x000000008000000aUL),
     SHA3_CONST(0x000000008000808bUL), SHA3_CONST(0x800000000000008bUL),
     SHA3_CONST(0x8000000000008089UL), SHA3_CONST(0x8000000000008003UL),
     SHA3_CONST(0x8000000000008002UL), SHA3_CONST(0x8000000000000080UL),
     SHA3_CONST(0x000000000000800aUL), SHA3_CONST(0x800000008000000aUL),
     SHA3_CONST(0x8000000080008081UL), SHA3_CONST(0x8000000000008080UL),
     SHA3_CONST(0x0000000080000001UL), SHA3_CONST(0x8000000080008008UL)
 };
 
 __device__ static const unsigned keccakf_rotc[24] = {
     1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14, 27, 41, 56, 8, 25, 43, 62,
     18, 39, 61, 20, 44
 };
 
 __device__ static const unsigned keccakf_piln[24] = {
     10, 7, 11, 17, 18, 3, 5, 16, 8, 21, 24, 4, 15, 23, 19, 13, 12, 2, 20,
     14, 22, 9, 6, 1
 };
 
 /* generally called after SHA3_KECCAK_SPONGE_WORDS-ctx->capacityWords words 
  * are XORed into the state s 
  */
  __device__ void keccakf(uint64_t s[25]) {
     int i, j, round;
     uint64_t t, bc[5];
     #define KECCAK_ROUNDS 24
 
     for (round = 0; round < KECCAK_ROUNDS; round++) {
         /* Theta */
         for (i = 0; i < 5; i++) {
             bc[i] = s[i] ^ s[i + 5] ^ s[i + 10] ^ s[i + 15] ^ s[i + 20];
         }
 
         for (i = 0; i < 5; i++) {
             t = bc[(i + 4) % 5] ^ SHA3_ROTL64(bc[(i + 1) % 5], 1);
             for (j = 0; j < 25; j += 5) {
                 s[j + i] ^= t;
             }
         }
 
         /* Rho Pi */
         t = s[1];
         for (i = 0; i < 24; i++) {
             j = keccakf_piln[i];
             bc[0] = s[j];
             s[j] = SHA3_ROTL64(t, keccakf_rotc[i]);
             t = bc[0];
         }
 
         /* Chi */
         for (j = 0; j < 25; j += 5) {
             for (i = 0; i < 5; i++) {
                 bc[i] = s[j + i];
             }
             for (i = 0; i < 5; i++) {
                 s[j + i] ^= (~bc[(i + 1) % 5]) & bc[(i + 2) % 5];
             }
         }
 
         /* Iota */
         s[0] ^= keccakf_rndc[round];
     }
 }
 