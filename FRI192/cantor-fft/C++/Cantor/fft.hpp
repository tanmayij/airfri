#ifndef ADDITIVE_FFT_CANTOR_HPP_
#define ADDITIVE_FFT_CANTOR_HPP_

#include <cstddef>
#include <vector>

#include "libiop/algebra/field_subset.hpp"
#include "libiop/algebra/subspace.hpp"

namespace cantor {

template<typename FieldT>
struct PreComputedValues {
    size_t dimension;
    std::vector<std::vector<size_t>> nz_S_vecs;
    std::vector<FieldT> mult_factor_vec;

    PreComputedValues(const size_t m){
        dimension = m;
        nz_S_vecs.reserve(m);
        mult_factor_vec.reserve((1<<m)-1);
    }
};

template<typename FieldT>
PreComputedValues<FieldT> pre_computation(const libiop::affine_subspace<FieldT> &domain);

template<typename FieldT>
std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                        const PreComputedValues<FieldT> &values);

template<typename FieldT>
std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &poly_coeffs,
                                        const PreComputedValues<FieldT> &values);

template<typename FieldT>
std::vector<FieldT> additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                        const libiop::affine_subspace<FieldT> &domain);

template<typename FieldT>
std::vector<FieldT> additive_IFFT(const std::vector<FieldT> &evals,
                                        const libiop::affine_subspace<FieldT> &domain);

} // namespace cantor

#include "Cantor/fft.tcc"
#endif // ADDITIVE_FFT_CANTOR_HPP_
