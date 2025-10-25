
//includes
#include <iostream>
#include <algorithm>
#include <set>
#include <vector>
#include <cstring>
#include <cstdint>
#include "../include/hash_host.cuh"
//sample 'number' unique indices in [0, size-1] for the first codeword, and fill reduced_indices as index % reduced_size
void sample_indices(uint8_t* seed, size_t seed_len, int size, int reduced_size, int number, size_t* indices, size_t* reduced_indices) {
    //zero-initialize output arrays
    std::fill(indices, indices + number, 0);
    std::fill(reduced_indices, reduced_indices + number, 0);
    std::cout << "[debug] size: " << size << ", number: " << number << std::endl;
    std::set<size_t> index_set;
    int counter = 0;
    uint8_t hash_output[32];
    while ((int)index_set.size() < number) {
        uint8_t input[seed_len + sizeof(int)];
        memcpy(input, seed, seed_len);
        memcpy(input + seed_len, &counter, sizeof(int));
        SHA3_host(hash_output, input, seed_len + sizeof(int), 256);
        uint64_t candidate = 0;
        for (int i = 0; i < 8; i++) {
            candidate = (candidate << 8) | hash_output[i];
        }
        size_t index = candidate % size;
        if (index_set.find(index) == index_set.end()) {
            index_set.insert(index);
        }
        counter++;
        if (counter > 1048577) {
            std::cerr << "error:could not sample enough unique indices" << std::endl;
            break;
        }
    }
    int idx = 0;
    for(auto it = index_set.begin(); it != index_set.end(); ++it) {
        indices[idx] = *it;
        reduced_indices[idx] = *it % reduced_size;
        idx++;
    }
    std::sort(indices, indices + number);
    std::sort(reduced_indices, reduced_indices + number);
    std::cout << "sampled initial indices: ";
    for(int i = 0; i < number; i++) {
        std::cout << indices[i] << " ";
    }
    std::cout << std::endl;
    std::cout << "reduced initial indices: ";
    for(int i = 0; i < number; i++) {
        std::cout << reduced_indices[i] << " ";
    }
    std::cout << std::endl;
}

//derive indices for the next round by dividing by 2 and deduplicating
void derive_next_round_indices(const size_t* prev_indices, int prev_number, size_t* next_indices, int& next_number) {
    std::set<size_t> next_set;
    for(int i = 0; i < prev_number; i++) {
        next_set.insert(prev_indices[i] / 2);
    }
    next_number = 0;
    for(auto it = next_set.begin(); it != next_set.end(); ++it) {
        next_indices[next_number++] = *it;
    }
    std::sort(next_indices, next_indices + next_number);
    std::cout << "next round indices: ";
    for(int i = 0; i < next_number; i++) {
        std::cout << next_indices[i] << " ";
    }
    std::cout << std::endl;
}