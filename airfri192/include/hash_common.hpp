#ifndef HASH_COMMON_HPP
#define HASH_COMMON_HPP

// Shared types, enums, and macros for SHA3/Keccak hash code (host and device)

#ifdef __cplusplus
extern "C" {
#endif

// SHA3 return codes
typedef enum {
    SHA3_OK = 0,
    SHA3_ERR = 1
} SHA3_RETURN;

// SHA3 flags
typedef enum {
    SHA3_FLAGS_NONE = 0,
    SHA3_FLAGS_KECCAK = 1
} SHA3_FLAGS;

// SHA3 context structure (minimal, extend as needed)
typedef struct {
    unsigned char saved[8];
    unsigned int s[50];
    unsigned int byteIndex;
    unsigned int wordIndex;
    unsigned int capacityWords;
    unsigned int fixedOutputLength;
    unsigned int bitsAvailableForSqueezing;
    int squeezing;
    SHA3_FLAGS flags;
} sha3_context;

#ifdef __cplusplus
}
#endif

#endif // HASH_COMMON_HPP
