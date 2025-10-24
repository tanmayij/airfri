#ifndef FIELD_CONVERSION_HPP
#define FIELD_CONVERSION_HPP

#include <cstdint>
#include <libff/algebra/fields/binary/gf256.hpp>
#include <libff/algebra/field_utils/field_utils.hpp>

// conversion between libff::gf256 and uint64_t[4] representation
// uint64_t[4] is used in CUDA kernels for compatibility with existing field operations

// convert libff::gf256 to uint64_t[4] (little-endian: word 0 = bits 0-63, word 3 = bits 192-255)
inline void field_to_uint64(const libff::gf256 &elem, uint64_t out[4]) {
    std::vector<uint64_t> words = elem.to_words();
    for (int i = 0; i < 4; ++i) out[i] = words[i];
}

// convert uint64_t[4] to libff::gf256 (little-endian)
inline libff::gf256 uint64_to_field(const uint64_t in[4]) {
    return libff::gf256(in[0], in[1], in[2], in[3]);
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
