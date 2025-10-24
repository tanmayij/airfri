#ifndef HASH_HPP
#define HASH_HPP

#include <cstdint>
#include <cstddef>

/* SHA-3/Keccak hash functions */

#define SHA3_KECCAK_SPONGE_WORDS \
    (((1600) / 8) / sizeof(uint64_t))

struct sha3_context {
    uint64_t saved;
    union {
        uint64_t s[SHA3_KECCAK_SPONGE_WORDS];
        uint8_t sb[SHA3_KECCAK_SPONGE_WORDS * 8];
    } u;
    unsigned byteIndex;
    unsigned wordIndex;
    unsigned capacityWords;
};

enum SHA3_FLAGS : uint32_t {
    SHA3_FLAGS_NONE = 0,
    SHA3_FLAGS_KECCAK = 1
};

enum SHA3_RETURN : uint32_t {
    SHA3_RETURN_OK = 0,
    SHA3_RETURN_BAD_PARAMS = 1
};

typedef enum SHA3_RETURN sha3_return_t;

/* Initialize SHA3 context */
sha3_return_t sha3_Init(void *priv, unsigned bitSize);
void sha3_Init256(void *priv);
void sha3_Init384(void *priv);
void sha3_Init512(void *priv);

// Inline host implementations
inline sha3_return_t sha3_Init(void *priv, unsigned bitSize) {
    sha3_context *ctx = (sha3_context *) priv;
    if( bitSize != 256 && bitSize != 384 && bitSize != 512 )
        return SHA3_RETURN_BAD_PARAMS;
    memset(ctx, 0, sizeof(*ctx));
    ctx->capacityWords = 2 * bitSize / (8 * sizeof(uint64_t));
    return SHA3_RETURN_OK;
}
inline void sha3_Init256(void *priv) { sha3_Init(priv, 256); }
inline void sha3_Init384(void *priv) { sha3_Init(priv, 384); }
inline void sha3_Init512(void *priv) { sha3_Init(priv, 512); }
inline enum SHA3_FLAGS sha3_SetFlags(void *priv, enum SHA3_FLAGS flags) {
    sha3_context *ctx = (sha3_context *) priv;
    flags = (enum SHA3_FLAGS)((uint32_t)flags & (uint32_t)SHA3_FLAGS_KECCAK);
    ctx->capacityWords |= (flags == SHA3_FLAGS_KECCAK ? 0x80000000 : 0);
    return flags;
}
inline void sha3_Update(void *priv, void const *bufIn, size_t len) {
    sha3_context *ctx = (sha3_context *) priv;
    unsigned old_tail = (8 - ctx->byteIndex) & 7;
    size_t words;
    unsigned tail;
    size_t i;
    const uint8_t *buf = (const uint8_t *)bufIn;
    if(len < old_tail) {
        while (len--)
            ctx->saved |= (uint64_t) (*(buf++)) << ((ctx->byteIndex++) * 8);
        return;
    }
    if(old_tail) {
        len -= old_tail;
        while (old_tail--)
            ctx->saved |= (uint64_t) (*(buf++)) << ((ctx->byteIndex++) * 8);
        ctx->u.s[ctx->wordIndex] ^= ctx->saved;
        ctx->byteIndex = 0;
        ctx->saved = 0;
        if(++ctx->wordIndex == (SHA3_KECCAK_SPONGE_WORDS - ctx->capacityWords)) {
            /* keccakf_host(ctx->u.s); */
            ctx->wordIndex = 0;
        }
    }
    words = len / sizeof(uint64_t);
    tail = len - words * sizeof(uint64_t);
    for(i = 0; i < words; i++, buf += sizeof(uint64_t)) {
        const uint64_t t = (uint64_t) (buf[0]) |
                ((uint64_t) (buf[1]) << 8 * 1) |
                ((uint64_t) (buf[2]) << 8 * 2) |
                ((uint64_t) (buf[3]) << 8 * 3) |
                ((uint64_t) (buf[4]) << 8 * 4) |
                ((uint64_t) (buf[5]) << 8 * 5) |
                ((uint64_t) (buf[6]) << 8 * 6) |
                ((uint64_t) (buf[7]) << 8 * 7);
        ctx->u.s[ctx->wordIndex] ^= t;
        if(++ctx->wordIndex == (SHA3_KECCAK_SPONGE_WORDS - ctx->capacityWords)) {
            /* keccakf_host(ctx->u.s); */
            ctx->wordIndex = 0;
        }
    }
    while (tail--) {
        ctx->saved |= (uint64_t) (*(buf++)) << ((ctx->byteIndex++) * 8);
    }
}
inline void const *sha3_Finalize(void *priv) {
    sha3_context *ctx = (sha3_context *) priv;
    uint64_t t;
    t = (uint64_t)(((uint64_t)(0x02 | (1 << 2))) << ((ctx->byteIndex) * 8));
    ctx->u.s[ctx->wordIndex] ^= ctx->saved ^ t;
    ctx->u.s[SHA3_KECCAK_SPONGE_WORDS - ctx->capacityWords - 1] ^= 0x8000000000000000UL;
    /* keccakf_host(ctx->u.s); */
    return (ctx->u.sb);
}
inline sha3_return_t sha3_HashBuffer(
    unsigned bitSize,
    enum SHA3_FLAGS flags,
    const void *in, unsigned inBytes,
    void *out, unsigned outBytes) {
    sha3_return_t err;
    sha3_context c;
    err = sha3_Init(&c, bitSize);
    if( err != SHA3_RETURN_OK )
        return err;
    if( sha3_SetFlags(&c, flags) != flags ) {
        return SHA3_RETURN_BAD_PARAMS;
    }
    sha3_Update(&c, in, inBytes);
    const void *h = sha3_Finalize(&c);
    if(outBytes > bitSize/8)
        outBytes = bitSize/8;
    memcpy(out, h, outBytes);
    return SHA3_RETURN_OK;
}
inline void SHA3s_host(uint8_t *hm, const uint8_t **msgs, size_t *msgs_len, size_t items, size_t bitSize)
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
inline void SHA3_host(uint8_t *hm, const uint8_t *msg, size_t msg_len, size_t bitSize)
{
    const uint8_t *msgs[1] = {msg};
    size_t lens[1] = {msg_len};
    SHA3s_host(hm, msgs, lens, 1, bitSize);
}
/* Set flags */
enum SHA3_FLAGS sha3_SetFlags(void *priv, enum SHA3_FLAGS flags);

/* Update with data */
void sha3_Update(void *priv, void const *bufIn, size_t len);

/* Finalize and get hash */
void const *sha3_Finalize(void *priv);

/* Single-call hash function */
void SHA3_host(uint8_t *hm, const uint8_t *msg, size_t msg_len, size_t bitSize);

/* Multiple messages hash */
void SHA3s_host(uint8_t *hm, const uint8_t **msgs, size_t *msgs_len, size_t items, size_t bitSize);

/* Buffer hash */
sha3_return_t sha3_HashBuffer(
    unsigned bitSize,
    enum SHA3_FLAGS flags,
    const void *in, unsigned inBytes,
    void *out, unsigned outBytes);

/* Convenience wrapper for SHA3-256 */
inline void hash_sha3_256(const uint8_t* input, size_t input_len, uint8_t* output) {
    SHA3_host(output, input, input_len, 256);
}

#endif // HASH_HPP
