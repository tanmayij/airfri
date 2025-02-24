// #include <stdio.h>
// #include <stdint.h>
// #include <string.h>
// #include "../include/hash-host.cuh"
// // Function prototype for your SHA3 implementation
// void SHA3_host(uint8_t *hm, const uint8_t *msg, size_t msg_len, size_t bitSize);

// int main() {
//     // Define the two uint64_t arrays
//     //89dd7af006338570 d48d9f6a0bd4789a 191a6d9b6064a1d9 99fcd542a536b8f0 e19f652cead15b7c 9b1d5b341353aefa 00d46fb70440a09a 673464b40434e931
//     uint64_t input1[] = {0x89dd7af006338570, 0xd48d9f6a0bd4789a, 0x191a6d9b6064a1d9, 0x99fcd542a536b8f0};
//     uint64_t input2[] = {0xe19f652cead15b7c, 0x9b1d5b341353aefa, 0x00d46fb70440a09a, 0x673464b40434e931};
//     //Computed Merkle Root: 7c77747851e67fa6 481cf7c1d470a73d 8e5c6e57f51975a5 f8b3b270e3457274 
//     // Convert input arrays to byte arrays for hashing
//     uint8_t msg[sizeof(input1) + sizeof(input2)];
//     memcpy(msg, input1, sizeof(input1));
//     memcpy(msg + sizeof(input1), input2, sizeof(input2));

//     // Output buffer for the SHA3-256 hash (32 bytes = 256 bits)
//     uint8_t hash[32];

//     // Compute SHA3-256 hash
//     SHA3_host(hash, msg, sizeof(msg), 256);

//     // Print the hash result
//     printf("SHA3-256 Hash: ");
//     for (size_t i = 0; i < sizeof(hash); i++) {
//         printf("%02x", hash[i]);
//     }
//     printf("\n");

//     return 0;
// }

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "../include/hash-host.cuh"  // Include your SHA3 implementation
#define HASH_WORDS 4
// Function to check if the given value produces the expected hash
int check_sha3_hash(uint64_t *input_value, uint64_t *expected_hash) {
    // Allocate memory for the computed hash
    uint64_t *computed_hash = (uint64_t *)calloc(HASH_WORDS,sizeof(uint64_t));
    if (computed_hash == NULL) {
        fprintf(stderr, "Memory allocation failed for computed_hash\n");
        exit(1);
    }

    // Allocate memory for intermediate hash storage
    uint64_t *new_hash = (uint64_t *)calloc(HASH_WORDS, sizeof(uint64_t));
    if (new_hash == NULL) {
        fprintf(stderr, "Memory allocation failed for new_hash\n");
        free(computed_hash);
        exit(1);
    }

    // Compute SHA3-256 hash
    printf("Input Bytes: ");
    uint8_t *byte_ptr = (uint8_t *)input_value;
    for (size_t i = 0; i < sizeof(uint64_t) * 2; i++) {
        printf("%02x ", byte_ptr[i]);
    }
    printf("\n");
    SHA3_host((uint8_t *)new_hash, (uint8_t *)input_value, 16 * sizeof(uint64_t), 256);

    // Copy new_hash into computed_hash
    memcpy(computed_hash, new_hash, HASH_WORDS * sizeof(uint64_t));

    // Reset new_hash to 0 for safety
    memset(new_hash, 0, HASH_WORDS * sizeof(uint64_t));

    // Debugging Output
    printf("Computed SHA3-256 Hash: ");
    for (int i = 0; i < HASH_WORDS; i++) {
        printf("%016lx ", computed_hash[i]);
    }
    printf("\n");

    // Compare computed hash with expected hash
    int result = memcmp(expected_hash, computed_hash, HASH_WORDS * sizeof(uint64_t)) == 0;

    // Free dynamically allocated memory
    free(computed_hash);
    free(new_hash);

    return result;
}

// Example usage
int main() {
    // Define sizes
    size_t input_size = 4 * 4 * sizeof(uint64_t);  // 4 uint64_t values (32 bytes)
    size_t hash_size = 4 * sizeof(uint64_t);   // SHA3-256 output (32 bytes)

    // Allocate memory for input values
    uint64_t *input_value = (uint64_t *)malloc(input_size);
    if (input_value == NULL) {
        fprintf(stderr, "Memory allocation failed for input_value\n");
        exit(1);
    }
    // uint64_t *input_value_two = (uint64_t *)malloc(input_size);
    // if (input_value_two == NULL) {
    //     fprintf(stderr, "Memory allocation failed for input_value\n");
    //     exit(1);
    // }
    // Assign values to input dynamically
    //aee4290157638f18 7376a1e1d60e9215 1ab57593ca9897c3 293169be78e5d870 52112ce2303df2ab 82e0aad5ebcb52dc 9d8330a46a595fe2 
    //300677d3a2bdd52b 8f4fcf245ff23125 f7946e66c2b269bb bbbeb802e252af6d 1aae94649c454267 c077c338c4e280ed 3fb2c191b54df591 3f77c65d96b80881 703ffcb6e4888675
    input_value[0] = 0xaee4290157638f18;
    input_value[1] = 0x7376a1e1d60e9215;
    input_value[2] = 0x1ab57593ca9897c3;
    input_value[3] = 0x293169be78e5d870;
    input_value[4] = 0x52112ce2303df2ab;
    input_value[5] = 0x82e0aad5ebcb52dc;
    input_value[6] = 0x9d8330a46a595fe2;
    input_value[7] = 0x300677d3a2bdd52b;
    input_value[8] = 0x8f4fcf245ff23125;
    input_value[9] = 0xf7946e66c2b269bb;
    input_value[10] = 0xbbbeb802e252af6d;
    input_value[11] = 0x1aae94649c454267;
    input_value[12] = 0xc077c338c4e280ed;
    input_value[13] = 0x3fb2c191b54df591;
    input_value[14] = 0x3f77c65d96b80881;
    input_value[15] = 0x703ffcb6e4888675;
   
    // Dynamically allocate memory for the expected hash
    uint64_t *expected_hash = (uint64_t *)malloc(hash_size);
    if (expected_hash == NULL) {
        fprintf(stderr, "Memory allocation failed for expected_hash\n");
        free(input_value);
        exit(1);
    }

    // Assign values to dynamically allocated hash array
    expected_hash[0] = 0x4ef977537b8a0509;
    expected_hash[1] = 0xede8eb9ea64a25c4;
    expected_hash[2] = 0x1cc96e1ab10ef543;
    expected_hash[3] = 0x839f06b4f98baba0;

    // Check if the input produces the expected hash
    int is_match = check_sha3_hash(input_value, expected_hash);
    //int is_is_match = check_sha3_hash(input_value_two, expected_hash);
    // Print the result
    if (is_match) {
        printf("Hash matches the expected output!\n");
    } else {
        printf("Hash does NOT match the expected output!\n");
    }

    // Free allocated memory
    free(input_value);
    free(expected_hash);
    //free(input_value_two);
    return 0;
}