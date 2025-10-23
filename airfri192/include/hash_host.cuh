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
