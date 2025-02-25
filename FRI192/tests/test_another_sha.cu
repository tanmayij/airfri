// #include <stdio.h>
// #include <stdint.h>
// #include <stdlib.h>
// #include <string.h>
// #include "../include/hash-host.cuh"  // Include your SHA3 implementation
// #define HASH_WORDS 4
// // Function to check if the given value produces the expected hash
// int check_sha3_hash(uint64_t *input_value, uint64_t *expected_hash) {
//     // Allocate memory for the computed hash
//     uint64_t *computed_hash = (uint64_t *)calloc(HASH_WORDS,sizeof(uint64_t));
//     if (computed_hash == NULL) {
//         fprintf(stderr, "Memory allocation failed for computed_hash\n");
//         exit(1);
//     }

//     // Allocate memory for intermediate hash storage
//     uint64_t *new_hash = (uint64_t *)calloc(HASH_WORDS, sizeof(uint64_t));
//     if (new_hash == NULL) {
//         fprintf(stderr, "Memory allocation failed for new_hash\n");
//         free(computed_hash);
//         exit(1);
//     }

//     // Compute SHA3-256 hash
//     printf("Input Bytes: ");
//     uint8_t *byte_ptr = (uint8_t *)input_value;
//     for (size_t i = 0; i < sizeof(uint64_t) * 2; i++) {
//         printf("%02x ", byte_ptr[i]);
//     }
//     printf("\n");
//     SHA3_host((uint8_t *)new_hash, (uint8_t *)input_value, 16 * sizeof(uint64_t), 256);

//     // Copy new_hash into computed_hash
//     memcpy(computed_hash, new_hash, HASH_WORDS * sizeof(uint64_t));

//     // Reset new_hash to 0 for safety
//     memset(new_hash, 0, HASH_WORDS * sizeof(uint64_t));

//     // Debugging Output
//     printf("Computed SHA3-256 Hash: ");
//     for (int i = 0; i < HASH_WORDS; i++) {
//         printf("%016lx ", computed_hash[i]);
//     }
//     printf("\n");

//     // Compare computed hash with expected hash
//     int result = memcmp(expected_hash, computed_hash, HASH_WORDS * sizeof(uint64_t)) == 0;

//     // Free dynamically allocated memory
//     free(computed_hash);
//     free(new_hash);

//     return result;
// }

// // Example usage
// int main() {
//     // Define sizes
//     size_t input_size = 4 * 2 * sizeof(uint64_t);  // 4 uint64_t values (32 bytes)
//     size_t hash_size = 4 * sizeof(uint64_t);   // SHA3-256 output (32 bytes)

//     // Allocate memory for input values
//     uint64_t *input_value = (uint64_t *)malloc(input_size);
//     if (input_value == NULL) {
//         fprintf(stderr, "Memory allocation failed for input_value\n");
//         exit(1);
//     }
//     // uint64_t *input_value_two = (uint64_t *)malloc(input_size);
//     // if (input_value_two == NULL) {
//     //     fprintf(stderr, "Memory allocation failed for input_value\n");
//     //     exit(1);
//     // }
//     // Assign values to input dynamically
//     //e5a7f5767e64120c 5489c9cff08d54af e1a9542d0fc7d723 909fe4aa0a71622b 29344544f00caa72 38f774874aeeb124 4cfa6381a6ad40b1 9ce4a8b4681cce0c
//     input_value[0] = 0xd46f02faa854a139;
//     input_value[1] = 0x5eda12fdcfc6b414;
//     input_value[2] = 0x4d2614da35978182;
//     input_value[3] = 0xd9c7afb2bf4dc4ae;
//     input_value[4] = 0xcc10bba2731e29a7;
//     input_value[5] = 0x4191eaf10b4fbbbf;
//     input_value[6] = 0xdd0931c636bd24bc;
//     input_value[7] = 0x35c56d2534530a85;
 
   
//     // Dynamically allocate memory for the expected hash
//     uint64_t *expected_hash = (uint64_t *)malloc(hash_size);
//     if (expected_hash == NULL) {
//         fprintf(stderr, "Memory allocation failed for expected_hash\n");
//         free(input_value);
//         exit(1);
//     }

//     // Assign values to dynamically allocated hash array
//     expected_hash[0] = 0xd46f02faa854a139;
//     expected_hash[1] = 0x5eda12fdcfc6b414;
//     expected_hash[2] = 0x4d2614da35978182;
//     expected_hash[3] = 0xd9c7afb2bf4dc4ae;

//     // Check if the input produces the expected hash
//     int is_match = check_sha3_hash(input_value, expected_hash);
//     //int is_is_match = check_sha3_hash(input_value_two, expected_hash);
//     // Print the result
//     if (is_match) {
//         printf("Hash matches the expected output!\n");
//     } else {
//         printf("Hash does NOT match the expected output!\n");
//     }

//     // Free allocated memory
//     free(input_value);
//     free(expected_hash);
//     //free(input_value_two);
//     return 0;
// }

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include "../include/hash-host.cuh"  // Ensure this includes the `SHA3_host` function

#define HASH_WORDS 4  // Assuming SHA3-256 output is 4 uint64_t values

void check_sha3_hash(const uint64_t input[HASH_WORDS], const uint64_t expected[HASH_WORDS]) {
    uint64_t output[HASH_WORDS];

    // Compute SHA3 hash
    SHA3_host((uint8_t *)output, (uint8_t *)input, HASH_WORDS * sizeof(uint64_t), 256);

    // Print computed hash
    printf("\nComputed SHA3-256 Hash:\n");
    for (int i = 0; i < HASH_WORDS; i++) {
        printf("%016lx ", output[i]);
    }
    printf("\n");

    // Print expected hash
    printf("Expected Hash:\n");
    for (int i = 0; i < HASH_WORDS; i++) {
        printf("%016lx ", expected[i]);
    }
    printf("\n");

    // Compare computed and expected values
    int match = 1;
    for (int i = 0; i < HASH_WORDS; i++) {
        if (output[i] != expected[i]) {
            match = 0;
            break;
        }
    }

    if (match) {
        printf("\n✅ Hash matches the expected output!\n");
    } else {
        printf("\n❌ Hash does NOT match the expected output!\n");
    }
}

// Example usage
int main() {
    uint64_t input[4 * HASH_WORDS] = {0x59adc3dd2ae6d02f, 0x7e7723533d8488ff, 0xe5fbf70703ff0535, 0xa238fef804d33cce, 0x44930ea537da00b2, 0x0e10126fe9483315, 0xb8c3b0f953cc0110, 0x8ecc941aedc4552d};
    uint64_t expected[HASH_WORDS] = {0xf3a531b17531e7ca, 0xc441865ce0816489, 0x88cd9af9aacef637, 0xee9149e44d8c2fd6};

    check_sha3_hash(input, expected);

    return 0;
}