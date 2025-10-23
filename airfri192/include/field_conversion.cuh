#ifndef FIELD_CONVERSION_HPP
#define FIELD_CONVERSION_HPP

#include <cstdint>
#include <libff/algebra/fields/binary/gf256.hpp>
#include <libff/algebra/field_utils/field_utils.hpp>

// conversion between libff::gf256 and uint64_t[4] representation
// uint64_t[4] is used in CUDA kernels for compatibility with existing field operations

// convert libff::gf256 to uint64_t[4] (little-endian: word 0 = bits 0-63, word 3 = bits 192-255)
inline void field_to_uint64(const libff::gf256 &elem, uint64_t out[4]) {
    libff::bit_vector bits = libff::convert_field_element_to_bit_vector(elem);
    
    for (int w = 0; w < 4; w++) {
        out[w] = 0;
        for (int b = 0; b < 64; b++) {
            if (bits[w * 64 + b]) {
                out[w] |= (1ULL << b);
            }
        }
    }
}

// convert uint64_t[4] to libff::gf256 (little-endian)
inline libff::gf256 uint64_to_field(const uint64_t in[4]) {
    libff::bit_vector bits(256);
    
    for (int w = 0; w < 4; w++) {
        for (int b = 0; b < 64; b++) {
            bits[w * 64 + b] = (in[w] >> b) & 1;
        }
    }
    
    return libff::convert_bit_vector_to_field_element<libff::gf256>(bits);
}

// batch conversion for vectors
inline void vector_to_uint64_array(const std::vector<libff::gf256> &vec, uint64_t *out) {
    for (size_t i = 0; i < vec.size(); i++) {
        field_to_uint64(vec[i], &out[i * 4]);
    }
}

// batch conversion from uint64_t array to vector
inline std::vector<libff::gf256> uint64_array_to_vector(const uint64_t *in, size_t num_elements) {
    std::vector<libff::gf256> result(num_elements);
    for (size_t i = 0; i < num_elements; i++) {
        result[i] = uint64_to_field(&in[i * 4]);
    }
    return result;
}

#endif // FIELD_CONVERSION_HPP
