#ifndef ADDITIVE_FFT_UTILS_HPP_
#define ADDITIVE_FFT_UTILS_HPP_

#include <cstddef>
#include <vector>


template<typename FieldT>
bool field_trace_binary(const FieldT &element);

template<typename FieldT>
std::vector<FieldT> cantor_basis(size_t m);

template <typename T>
void my_print_vector(const std::vector<T> &v);

template <typename FieldT>
FieldT divide(std::vector<FieldT> &g, const std::vector<size_t> &nz_S, const size_t input_size, const size_t offset); 

#include "utils/utils.tcc"

#endif // ADDITIVE_FFT_UTILS_HPP_
