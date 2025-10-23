#include "sample_indices.hpp"
#include "hash.hpp"
#include <iostream>
#include <cstring>
#include <algorithm>

void sample_indices(
    uint8_t* seed, 
    size_t seed_len, 
    int size, 
    int reduced_size, 
    int number, 
    size_t* indices, 
    size_t* reduced_indices
) {
    // Validate inputs
    if (number > reduced_size) {
        std::cerr << "Error: number > reduced_size" << std::endl;
        return;
    }
    if (number > 2 * reduced_size) {
        std::cerr << "Error: number > 2 * reduced_size" << std::endl;
        return;
    }
    
    std::cout << "Sampling indices with:" << std::endl;
    std::cout << "  Seed (first 8 bytes): ";
    for (size_t i = 0; i < std::min(seed_len, (size_t)8); i++) {
        printf("%02x ", seed[i]);
    }
    std::cout << std::endl;
    std::cout << "  size: " << size << std::endl;
    std::cout << "  reduced_size: " << reduced_size << std::endl;
    std::cout << "  number: " << number << std::endl;
    
    int counter = 0;
    int indices_count = 0;
    int reduced_indices_count = 0;
    
    uint8_t hash_output[32];
    
    while (indices_count < number) {
        // Hash seed || counter to get deterministic randomness
        uint8_t input[seed_len + sizeof(int)];
        memcpy(input, seed, seed_len);
        memcpy(input + seed_len, &counter, sizeof(int));
        
        // Hash to get 32 bytes of pseudorandom data
        SHA3_host(hash_output, input, seed_len + sizeof(int), 256);
        
        // Extract 4-byte chunks as potential indices
        for (int i = 0; i < 32; i += 4) {
            if (indices_count >= number) break;
            
            // Interpret 4 bytes as uint32_t
            uint32_t candidate = 0;
            for (int j = 0; j < 4; j++) {
                candidate |= ((uint32_t)hash_output[i + j]) << (j * 8);
            }
            
            // Map to range [0, size)
            size_t index = candidate % size;
            
            // Check for duplicates
            bool duplicate = false;
            for (int k = 0; k < indices_count; k++) {
                if (indices[k] == index) {
                    duplicate = true;
                    break;
                }
            }
            
            if (!duplicate) {
                indices[indices_count] = index;
                
                // Compute reduced index (for next round)
                size_t reduced_index = index % reduced_size;
                
                // Check if reduced index is duplicate
                bool reduced_duplicate = false;
                for (int k = 0; k < reduced_indices_count; k++) {
                    if (reduced_indices[k] == reduced_index) {
                        reduced_duplicate = true;
                        break;
                    }
                }
                
                if (!reduced_duplicate && reduced_indices_count < number) {
                    reduced_indices[reduced_indices_count] = reduced_index;
                    reduced_indices_count++;
                }
                
                indices_count++;
            }
        }
        
        counter++;
        
        // Safety check to prevent infinite loop
        if (counter > 10000) {
            std::cerr << "Error: Could not sample enough unique indices" << std::endl;
            break;
        }
    }
    
    // Sort indices for consistent ordering
    std::sort(indices, indices + indices_count);
    std::sort(reduced_indices, reduced_indices + reduced_indices_count);
    
    // Print sampled indices
    std::cout << "Sampled indices: ";
    for (int i = 0; i < indices_count; i++) {
        std::cout << indices[i] << " ";
    }
    std::cout << std::endl;
    
    std::cout << "Reduced indices: ";
    for (int i = 0; i < reduced_indices_count; i++) {
        std::cout << reduced_indices[i] << " ";
    }
    std::cout << std::endl;
}
