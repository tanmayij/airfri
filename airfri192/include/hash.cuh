#ifndef SHA3_H
#define SHA3_H

#include <cstdint>
#include <cstddef>

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

/* 'Words' here refers to uint64_t */
#define SHA3_KECCAK_SPONGE_WORDS \
    (((1600) / 8 /*bits to byte*/) / sizeof(uint64_t))
typedef struct sha3_context_
{
    uint64_t saved; /* the portion of the input message that we
                     * didn't consume yet */
    union
    { /* Keccak's state */
        uint64_t s[SHA3_KECCAK_SPONGE_WORDS];
        uint8_t sb[SHA3_KECCAK_SPONGE_WORDS * 8];
    } u;
    unsigned byteIndex;     /* 0..7--the next byte after the set one
                             * (starts from 0; 0--none are buffered) */
    unsigned wordIndex;     /* 0..24--the next word to integrate input
                             * (starts from 0) */
    unsigned capacityWords; /* the double size of the hash output in
                             * words (e.g. 16 for Keccak 512) */
} sha3_context;

enum SHA3_FLAGS : uint32_t
{
    SHA3_FLAGS_NONE = 0,
    SHA3_FLAGS_KECCAK = 1
};

enum SHA3_RETURN: uint32_t
{
    SHA3_RETURN_OK = 0,
    SHA3_RETURN_BAD_PARAMS = 1
};
typedef enum SHA3_RETURN sha3_return_t;

/* For Init or Reset call these: */

#ifdef __CUDACC__
// Inline device implementations
inline __device__ sha3_return_t sha3_Init(void *priv, unsigned bitSize) {
    sha3_context *ctx = (sha3_context *)priv;
    if (bitSize != 256 && bitSize != 384 && bitSize != 512) {
        return SHA3_RETURN_BAD_PARAMS;
    }
    memset(ctx, 0, sizeof(*ctx));
    ctx->capacityWords = 2 * bitSize / (8 * sizeof(uint64_t));
    return SHA3_RETURN_OK;
}

inline __device__ void sha3_Init256(void *priv) { sha3_Init(priv, 256); }
inline __device__ void sha3_Init384(void *priv) { sha3_Init(priv, 384); }
inline __device__ void sha3_Init512(void *priv) { sha3_Init(priv, 512); }

inline __device__ enum SHA3_FLAGS sha3_SetFlags(void *priv, enum SHA3_FLAGS flags) {
    sha3_context *ctx = (sha3_context *)priv;
    flags = (enum SHA3_FLAGS)((uint32_t)flags & (uint32_t)SHA3_FLAGS_KECCAK);
    ctx->capacityWords |= (flags == SHA3_FLAGS_KECCAK ? 0x80000000 : 0);
    return flags;
}

inline __device__ void sha3_Update(void *priv, void const *bufIn, size_t len) {
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
        if (++ctx->wordIndex == (SHA3_KECCAK_SPONGE_WORDS - ctx->capacityWords)) {
            /* keccakf(ctx->u.s); */
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
        if (++ctx->wordIndex == (SHA3_KECCAK_SPONGE_WORDS - ctx->capacityWords)) {
            /* keccakf(ctx->u.s); */
            ctx->wordIndex = 0;
        }
    }
    while (tail--) {
        ctx->saved |= (uint64_t)(*(buf++)) << ((ctx->byteIndex++) * 8);
    }
}

inline __device__ void const *sha3_Finalize(void *priv) {
    sha3_context *ctx = (sha3_context *)priv;
    uint64_t t;
    t = (uint64_t)(((uint64_t)(0x02 | (1 << 2))) << ((ctx->byteIndex) * 8));
    ctx->u.s[ctx->wordIndex] ^= ctx->saved ^ t;
    ctx->u.s[SHA3_KECCAK_SPONGE_WORDS - ctx->capacityWords - 1] ^= 0x8000000000000000UL;
    /* keccakf(ctx->u.s); */
    return (ctx->u.sb);
}

inline __device__ sha3_return_t sha3_HashBuffer(
    unsigned bitSize,
    enum SHA3_FLAGS flags,
    const void *in, unsigned inBytes,
    void *out, unsigned outBytes) {
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

inline __device__ void SHA3s(uint8_t *hm, const uint8_t **msgs, size_t *msgs_len, size_t items, size_t bitSize)
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
inline __device__ void SHA3(uint8_t *hm, const uint8_t *msg, size_t msg_len, size_t bitSize)
{
    SHA3s(hm, (const uint8_t *[]){msg}, (size_t[]){msg_len}, 1, bitSize);
}
#endif

#endif