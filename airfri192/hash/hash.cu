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
 #include <stdio.h>
 #include <stdint.h>
 #include <string.h>
 
 #include "hash.cuh"
 
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
 
 /* *************************** Public Interface ************************ */
 
 /* For Init or Reset call these: */
 __device__ sha3_return_t sha3_Init(void *priv, unsigned bitSize) {
     sha3_context *ctx = (sha3_context *)priv;
     if (bitSize != 256 && bitSize != 384 && bitSize != 512) {
         return SHA3_RETURN_BAD_PARAMS;
     }
     memset(ctx, 0, sizeof(*ctx));  // Replace with cuda-compatible memset if needed
     ctx->capacityWords = 2 * bitSize / (8 * sizeof(uint64_t));
     return SHA3_RETURN_OK;
 }
 
 __device__ void sha3_Init256(void *priv) {
     sha3_Init(priv, 256);
 }
 
 __device__ void sha3_Init384(void *priv) {
     sha3_Init(priv, 384);
 }
 
 __device__ void sha3_Init512(void *priv) {
     sha3_Init(priv, 512);
 }
 
 __device__ enum SHA3_FLAGS sha3_SetFlags(void *priv, enum SHA3_FLAGS flags) {
     sha3_context *ctx = (sha3_context *)priv;
     flags = (enum SHA3_FLAGS)((uint32_t)flags & (uint32_t)SHA3_FLAGS_KECCAK);
     ctx->capacityWords |= (flags == SHA3_FLAGS_KECCAK ? SHA3_USE_KECCAK_FLAG : 0);
     return flags;
 }
 
 __device__ void sha3_Update(void *priv, void const *bufIn, size_t len) {
     sha3_context *ctx = (sha3_context *)priv;
 
     unsigned old_tail = (8 - ctx->byteIndex) & 7;
 
     size_t words;
     unsigned tail;
     size_t i;
 
     const uint8_t *buf = (const uint8_t *)bufIn;
 
     if (len < old_tail) {
         while (len--) {
             ctx->saved |= (uint64_t)(*(buf++)) << ((ctx->byteIndex++) * 8);
         }
         return;
     }
 
     if (old_tail) {
         len -= old_tail;
         while (old_tail--) {
             ctx->saved |= (uint64_t)(*(buf++)) << ((ctx->byteIndex++) * 8);
         }
         ctx->u.s[ctx->wordIndex] ^= ctx->saved;
         ctx->byteIndex = 0;
         ctx->saved = 0;
         if (++ctx->wordIndex == (SHA3_KECCAK_SPONGE_WORDS - SHA3_CW(ctx->capacityWords))) {
             keccakf(ctx->u.s);
             ctx->wordIndex = 0;
         }
     }
 
     words = len / sizeof(uint64_t);
     tail = len - words * sizeof(uint64_t);
 
     for (i = 0; i < words; i++, buf += sizeof(uint64_t)) {
         uint64_t t = (uint64_t)(buf[0]) |
                      ((uint64_t)(buf[1]) << 8 * 1) |
                      ((uint64_t)(buf[2]) << 8 * 2) |
                      ((uint64_t)(buf[3]) << 8 * 3) |
                      ((uint64_t)(buf[4]) << 8 * 4) |
                      ((uint64_t)(buf[5]) << 8 * 5) |
                      ((uint64_t)(buf[6]) << 8 * 6) |
                      ((uint64_t)(buf[7]) << 8 * 7);
         ctx->u.s[ctx->wordIndex] ^= t;
         if (++ctx->wordIndex == (SHA3_KECCAK_SPONGE_WORDS - SHA3_CW(ctx->capacityWords))) {
             keccakf(ctx->u.s);
             ctx->wordIndex = 0;
         }
     }
 
     while (tail--) {
         ctx->saved |= (uint64_t)(*(buf++)) << ((ctx->byteIndex++) * 8);
     }
 }
 
 __device__ void const *sha3_Finalize(void *priv) {
     sha3_context *ctx = (sha3_context *)priv;
     uint64_t t;
 
     if (ctx->capacityWords & SHA3_USE_KECCAK_FLAG) {
         t = (uint64_t)(((uint64_t)1) << (ctx->byteIndex * 8));
     } else {
         t = (uint64_t)(((uint64_t)(0x02 | (1 << 2))) << ((ctx->byteIndex) * 8));
     }
 
     ctx->u.s[ctx->wordIndex] ^= ctx->saved ^ t;
     ctx->u.s[SHA3_KECCAK_SPONGE_WORDS - SHA3_CW(ctx->capacityWords) - 1] ^= SHA3_CONST(0x8000000000000000UL);
     keccakf(ctx->u.s);
 
     return (ctx->u.sb);
 }
 
 __device__ sha3_return_t sha3_HashBuffer(unsigned bitSize, enum SHA3_FLAGS flags, const void *in, unsigned inBytes, void *out, unsigned outBytes) {
     sha3_return_t err;
     sha3_context c;
 
     err = sha3_Init(&c, bitSize);
     if (err != SHA3_RETURN_OK) {
         return err;
     }
     if (sha3_SetFlags(&c, flags) != flags) {
         return SHA3_RETURN_BAD_PARAMS;
     }
     sha3_Update(&c, in, inBytes);
     const void *h = sha3_Finalize(&c);
 
     if (outBytes > bitSize / 8) {
         outBytes = bitSize / 8;
     }
     memcpy(out, h, outBytes);
     return SHA3_RETURN_OK;
 }

extern "C" __device__ void SHA3(uint8_t *hm, const uint8_t *msg, size_t msg_len, size_t bitSize)
{   
    SHA3s(hm, (const uint8_t *[]){msg}, (size_t[]){msg_len}, 1, bitSize);
}

__device__ void SHA3s(uint8_t *hm, const uint8_t **msgs, size_t *msgs_len, size_t items, size_t bitSize)
{
    sha3_context ctx;
    sha3_Init(&ctx, bitSize);
    for (size_t i = 0; i < items; ++i)
    {
        sha3_Update(&ctx, msgs[i], msgs_len[i]);
    }
    const void *result = sha3_Finalize(&ctx);
    memcpy(hm, result, bitSize / 8);
}