#include <iostream>
#include <vector>
#include <cstring>
#include <cstdlib>
#include <openssl/evp.h>

// Define a maximum number of objects for safety
#define MAX_PROOF_STREAM_OBJECTS 131072

unsigned char *global_serialized_data = nullptr;
size_t global_serialized_size = 0;
size_t global_prover_serialized_position = 0;
size_t global_serialized_position = 0;

// Define a structure for the ProofStream
struct ProofStream {
    std::vector<void*> objects;
    size_t read_index;

    ProofStream() : read_index(0) {
        objects.reserve(1024); // Reserve some initial capacity
    }

    ~ProofStream() {
        for (void* obj : objects) {
            free(obj);
        }
    }
};

// Initialize a ProofStream
ProofStream* create_proof_stream() {
    return new ProofStream();
}

// Push an object into the ProofStream with a maximum limit check
void push_object(ProofStream *ps, void *obj) {
    if (ps->objects.size() >= MAX_PROOF_STREAM_OBJECTS) {
        std::cerr << "Error: ProofStream has reached its maximum capacity of " 
                  << MAX_PROOF_STREAM_OBJECTS << " objects." << std::endl;
        exit(EXIT_FAILURE);
    }

    ps->objects.push_back(obj);
}

// Pull an object from the ProofStream
void* pull_object(ProofStream *ps) {
    if (ps->read_index < ps->objects.size()) {
        return ps->objects[ps->read_index++];
    } else {
        std::cerr << "ProofStream: cannot pull object; queue empty." << std::endl;
        exit(EXIT_FAILURE);
    }
}

// Serialize the ProofStream
unsigned char* serialize_proof_stream(ProofStream *ps, size_t num_bytes, size_t *out_size) {
    const size_t HASH_SIZE = num_bytes;
    size_t total_size = ps->objects.size() * HASH_SIZE;
    
    unsigned char *buffer = (unsigned char*) malloc(total_size);
    if (buffer == nullptr) {
        std::cerr << "Failed to allocate memory for serialization" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    // Copy each object into the buffer
    unsigned char *current_position = buffer;
    for (void* obj : ps->objects) {
        memcpy(current_position, obj, HASH_SIZE);
        current_position += HASH_SIZE;
    }
    
    *out_size = total_size;
    return buffer;
}

// Deserialize the ProofStream
ProofStream* deserialize_proof_stream(unsigned char *data, size_t size) {
    ProofStream *ps = create_proof_stream();
    size_t num_objects = size / sizeof(void*);
    
    for (size_t i = 0; i < num_objects; i++) {
        void* obj = malloc(sizeof(void*));
        if (obj == nullptr) {
            std::cerr << "Failed to allocate memory for deserialization" << std::endl;
            exit(EXIT_FAILURE);
        }
        memcpy(obj, data + i * sizeof(void*), sizeof(void*));
        ps->objects.push_back(obj);
    }
    
    ps->read_index = 0;
    return ps;
}

// Helper function to compute SHAKE-256
unsigned char* compute_shake256(unsigned char *input, size_t input_len, size_t num_bytes) {
    EVP_MD_CTX *mdctx;
    unsigned char *md_value = (unsigned char*) malloc(num_bytes);
    if (md_value == nullptr) {
        std::cerr << "Failed to allocate memory for hash" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    mdctx = EVP_MD_CTX_new();
    EVP_DigestInit_ex(mdctx, EVP_shake256(), nullptr);
    EVP_DigestUpdate(mdctx, input, input_len);
    EVP_DigestFinalXOF(mdctx, md_value, num_bytes);
    EVP_MD_CTX_free(mdctx);
    
    return md_value;
}

// Prover's Fiat-Shamir transformation
unsigned char* prover_fiat_shamir(ProofStream *ps, size_t num_bytes) {
    size_t serialized_size = 0;
    unsigned char *serialized_data = serialize_proof_stream(ps, num_bytes, &serialized_size);
    
    unsigned char *hash = compute_shake256(
        serialized_data + global_prover_serialized_position, 
        num_bytes, 
        num_bytes
    );
    
    global_prover_serialized_position += num_bytes;
    
    if (global_serialized_data != nullptr) {
        free(global_serialized_data);
    }
    global_serialized_data = serialized_data;
    global_serialized_size = serialized_size;
    
    return hash;
}

// Verifier's Fiat-Shamir transformation
unsigned char* verifier_fiat_shamir(ProofStream *ps, size_t num_bytes) {
    if (global_serialized_data == nullptr) {
        std::cerr << "Global serialized data is not set by the prover." << std::endl;
        exit(EXIT_FAILURE);
    }
    
    if (global_serialized_position + num_bytes > global_serialized_size) {
        std::cerr << "Not enough data in the global serialized buffer." << std::endl;
        exit(EXIT_FAILURE);
    }

    unsigned char *temp_buffer = (unsigned char*) malloc(num_bytes);
    if (temp_buffer == nullptr) {
        std::cerr << "Failed to allocate memory for temp buffer" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    memcpy(temp_buffer, global_serialized_data + global_serialized_position, num_bytes);
    global_serialized_position += num_bytes;

    unsigned char *hash = compute_shake256(temp_buffer, num_bytes, num_bytes);
    free(temp_buffer);
    
    return hash;
}

// Free the proof stream
void free_proof_stream(ProofStream* ps) {
    delete ps;
}

// Clear global serialized data
void clear_global_serialized_data() {
    if (global_serialized_data != nullptr) {
        free(global_serialized_data);
        global_serialized_data = nullptr;
    }
    global_serialized_size = 0;
    global_serialized_position = 0;
    global_prover_serialized_position = 0;
}
