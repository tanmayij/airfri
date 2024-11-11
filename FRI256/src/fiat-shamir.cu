#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <openssl/evp.h>

// Define a maximum number of objects for safety
#define MAX_PROOF_STREAM_OBJECTS 2097152

unsigned char *global_serialized_data = NULL;
size_t global_serialized_size = 0;
size_t global_prover_serialized_position = 0;
size_t global_serialized_position = 0;

// Define a structure for the ProofStream
typedef struct {
    void **objects;
    size_t num_objects;
    size_t read_index;
} ProofStream;

// Initialize a ProofStream
ProofStream* create_proof_stream() {
    ProofStream *ps = (ProofStream*) malloc(sizeof(ProofStream));
    ps->objects = (void **)malloc(sizeof(void *));
    ps->num_objects = 0;
    ps->read_index = 0;
    return ps;
}

// Push an object into the ProofStream with a maximum limit check
void push_object(ProofStream *ps, void *obj) {
    if (ps->num_objects >= MAX_PROOF_STREAM_OBJECTS) {
        fprintf(stderr, "Error: ProofStream has reached its maximum capacity of %d objects.\n", MAX_PROOF_STREAM_OBJECTS);
        exit(EXIT_FAILURE); // Exit or handle error as desired
    }

    ps->objects = (void**)realloc(ps->objects, sizeof(void*) * (ps->num_objects + 1));
    if (ps->objects == NULL) {
        fprintf(stderr, "Error: Failed to allocate memory for ProofStream objects.\n");
        exit(EXIT_FAILURE); // Exit or handle error as desired
    }

    ps->objects[ps->num_objects] = obj;
    ps->num_objects++;
    //printf("ProofStream after push: num_objects = %zu, read_index = %zu\n", ps->num_objects, ps->read_index);
}
// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <openssl/evp.h>

// unsigned char *global_serialized_data = NULL;
// size_t global_serialized_size = 0;
// size_t global_prover_serialized_position = 0;
// size_t global_serialized_position = 0;
// // Define a structure for the ProofStream
// typedef struct {
//     void **objects;
//     size_t num_objects;
//     size_t read_index;
// } ProofStream;

// // Initialize a ProofStream
// ProofStream* create_proof_stream() {
//     ProofStream *ps = (ProofStream*) malloc(sizeof(ProofStream));
//     ps->objects = NULL;
//     ps->num_objects = 0;
//     ps->read_index = 0;
//     return ps;
// }

// // Push an object into the ProofStream
// void push_object(ProofStream *ps, void *obj) {
//     ps->objects = (void**)realloc(ps->objects, sizeof(void*) * (ps->num_objects + 1));
//     ps->objects[ps->num_objects] = obj;
//     ps->num_objects++;
//     printf("ProofStream after push: num_objects = %zu, read_index = %zu\n", ps->num_objects, ps->read_index);
//     // for (size_t i = 0; i < ps->num_objects; i++) {
//     //     printf("  Object %zu: %p\n", i, ps->objects[ps->num_objects++]);
//     // }
// }

// Pull an object from the ProofStream
void* pull_object(ProofStream *ps) {
    if (ps->read_index < ps->num_objects) {
        //printf("ProofStream after pull: num_objects = %zu, read_index = %zu\n", ps->num_objects, ps->read_index);
        // for (size_t i = 0; i < ps->num_objects; i++) {
        // printf("  Object %zu: %p\n", i, ps->objects[ps->num_objects++]);
        // }
        return ps->objects[ps->read_index++];
        //printf("inside fiat-shamir.c: %p ", ps->objects[ps->read_index++]);
    } else {
        fprintf(stderr, "ProofStream: cannot pull object; queue empty.\n");
        exit(EXIT_FAILURE);
    }
}

// // Helper function to serialize the ProofStream (basic version)
// unsigned char* serialize_proof_stream(ProofStream *ps, size_t *out_size) {
//     size_t total_size = ps->num_objects * sizeof(void*);
//     unsigned char *buffer = (unsigned char*) malloc(total_size);
//     if (buffer == NULL) {
//         fprintf(stderr, "Failed to allocate memory for serialization\n");
//         exit(EXIT_FAILURE);
//     }
//     memcpy(buffer, ps->objects, total_size);
//     *out_size = total_size;
//     return buffer;
// }
unsigned char* serialize_proof_stream(ProofStream *ps, size_t num_bytes, size_t *out_size) {
    const size_t HASH_SIZE = num_bytes;

    size_t total_size = ps->num_objects * HASH_SIZE;
    //printf("total_size: %zu \n", total_size);
    unsigned char *buffer = (unsigned char*) malloc(total_size);
    if (buffer == NULL) {
        fprintf(stderr, "Failed to allocate memory for serialization\n");
        exit(EXIT_FAILURE);
    }
    //printf("number of objects in fiat-shamir serialization? %zu\n", ps->num_objects);
    // Copy each object into the buffer
    unsigned char *current_position = buffer;
    for (size_t i = 0; i < ps->num_objects; i++) {
        memcpy(current_position, ps->objects[i], HASH_SIZE);
        current_position += HASH_SIZE;
    }
    *out_size = total_size;
    // for (size_t i = 0; i < total_size; i++) {
    //     printf("%02x ", buffer[i]);
    // }
    // printf("buffer\n");
    return buffer;
}

ProofStream* deserialize_proof_stream(unsigned char *data, size_t size) {
    ProofStream *ps = create_proof_stream();
    size_t num_objects = size / sizeof(void*);
    ps->objects = (void**) malloc(num_objects * sizeof(void*));
    if (ps->objects == NULL) {
        fprintf(stderr, "Failed to allocate memory for deserialization\n");
        exit(EXIT_FAILURE);
    }
    memcpy(ps->objects, data, size);
    ps->num_objects = num_objects;
    ps->read_index = 0;
    return ps;
}

// Helper function to compute SHAKE-256
unsigned char* compute_shake256(unsigned char *input, size_t input_len, size_t num_bytes) {
    EVP_MD_CTX *mdctx;
    unsigned char *md_value = (unsigned char*) malloc(num_bytes);
    if (md_value == NULL) {
        fprintf(stderr, "Failed to allocate memory for hash\n");
        exit(EXIT_FAILURE);
    }
    mdctx = EVP_MD_CTX_new();
    EVP_DigestInit_ex(mdctx, EVP_shake256(), NULL);
    EVP_DigestUpdate(mdctx, input, input_len);
    EVP_DigestFinalXOF(mdctx, md_value, num_bytes);
    EVP_MD_CTX_free(mdctx);
    return md_value;
}

// Prover's Fiat-Shamir transformation
unsigned char* prover_fiat_shamir(ProofStream *ps, size_t num_bytes) {
    size_t serialized_size = ps->read_index * sizeof(void*); //this just sets to 0, but there seems to be no issue with that
    //printf("serialized size: %zu \n", serialized_size);
    unsigned char *serialized_data = serialize_proof_stream(ps, num_bytes, &serialized_size);
    // Print serialized_data
    // printf("Prover Serialized data after memcpy: ");
    // for (size_t i = 0; i < serialized_size; i++) {
    //     printf("%02x ", serialized_data[i]);
    // }
    // printf("\n");
    unsigned char *hash = compute_shake256(serialized_data + global_prover_serialized_position, num_bytes, num_bytes);
    global_prover_serialized_position += num_bytes;
    global_serialized_data = serialized_data;
    global_serialized_size = serialized_size;
    free(serialized_data);
    return hash;
}


//Verifier's Fiat-Shamir transformation
// unsigned char* verifier_fiat_shamir(ProofStream *ps, size_t num_bytes) {
//     size_t serialized_size = ps->read_index * sizeof(void*);
//     unsigned char *serialized_data = (unsigned char*) malloc(serialized_size);
//     if (serialized_data == NULL) {
//         fprintf(stderr, "Failed to allocate memory for serialization\n");
//         exit(EXIT_FAILURE);
//     }
//     // Print objects before memcpy
//     printf("ProofStream objects before memcpy:\n");
//     for (size_t i = 0; i < ps->read_index; i++) {
//         uint64_t *obj = (uint64_t *)ps->objects[i];
//         printf("Object %zu: ", i);
//         for (size_t j = 0; j < 4; j++) { // Assuming each object is 4 uint64_t
//             printf("%016llx ", obj[j]);
//         }
//         printf("\n");
//     }

//     memcpy(serialized_data, ps->objects, serialized_size);

//     unsigned char *hash = compute_shake256(serialized_data, serialized_size, num_bytes);
//     free(serialized_data);
//     return hash;
// }
unsigned char* verifier_fiat_shamir(ProofStream *ps, size_t num_bytes) {
    if (global_serialized_data == NULL) {
        fprintf(stderr, "Global serialized data is not set by the prover.\n");
        exit(EXIT_FAILURE);
    }
    if (global_serialized_position + num_bytes > global_serialized_size) {
        fprintf(stderr, "Not enough data in the global serialized buffer.\n");
        exit(EXIT_FAILURE);
    }

    unsigned char *temp_buffer = (unsigned char*) malloc(num_bytes);
    if (temp_buffer == NULL) {
        fprintf(stderr, "Failed to allocate memory for temp buffer\n");
        exit(EXIT_FAILURE);
    }
    
    memcpy(temp_buffer, global_serialized_data + global_serialized_position, num_bytes);
    global_serialized_position += num_bytes;

    // Print temp_buffer for debugging
    // printf("Verifier Serialized data: ");
    // for (size_t i = 0; i < num_bytes; i++) {
    //     printf("%02x ", temp_buffer[i]);
    // }
    // printf("\n");

    unsigned char *hash = compute_shake256(temp_buffer, num_bytes, num_bytes);
    free(temp_buffer);
    return hash;
}


void free_proof_stream(ProofStream* ps) {
    for (size_t i = 0; i < ps->num_objects; i++) {
        free(ps->objects[i]);
    }
    free(ps->objects);
    free(ps);
}
void clear_global_serialized_data() {
    if (global_serialized_data != NULL) {
        free(global_serialized_data);
        global_serialized_data = NULL;
    }
    global_serialized_size = 0;
    global_serialized_position = 0;

    printf("Global serialized data cleared.\n");
}
// int main() {
//     ProofStream *ps = create_proof_stream();
//     int obj1 = 42; // Example object
//     int obj2 = 24; // Another example object
//     int obj3 = 64; // Another example object
//     push_object(ps, &obj1);
//     push_object(ps, &obj2);
//     push_object(ps, &obj3);

//     // Debugging: Print proof_stream details
//     printf("Main Function: proof_stream->num_objects = %zu\n", ps->num_objects);
//     printf("Main Function: proof_stream->read_index = %zu\n", ps->read_index);

//     unsigned char *prover_hash = prover_fiat_shamir(ps, 32);
//     unsigned char *verifier_hash = verifier_fiat_shamir(ps, 32);

//     printf("Prover Hash: ");
//     for (int i = 0; i < 32; i++) printf("%02x", prover_hash[i]);
//     printf("\n");

//     printf("Verifier Hash: ");
//     for (int i = 0; i < 32; i++) printf("%02x", verifier_hash[i]);
//     printf("\n");

//     printf("Pulling objects from ProofStream:\n");
//     int *pulled_obj1 = (int*)pull_object(ps);
//     int *pulled_obj2 = (int*)pull_object(ps);
//     int *pulled_obj3 = (int*)pull_object(ps);
//     printf("Pulled Object 1: %d\n", *pulled_obj1);
//     printf("Pulled Object 2: %d\n", *pulled_obj2);
//     printf("Pulled Object 3: %d\n", *pulled_obj3);

//     free(prover_hash);
//     free(verifier_hash);
//     free_proof_stream(ps);

//     return 0;
// }