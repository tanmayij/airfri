#ifndef ADDITIVE_FFT_GAO_HPP_
#define ADDITIVE_FFT_GAO_HPP_

#include <cstddef>
#include <vector>

#include "libiop/algebra/field_subset.hpp"
#include "libiop/algebra/subspace.hpp"

namespace gao {

template <typename FieldT> 
struct PreComputedValues_Level1 {

  std::vector<FieldT> mult_factor_vec;
  std::vector<FieldT> recursed_betas;
  std::vector<FieldT> recursed_shifts;
  size_t dimension;

  PreComputedValues_Level1(const size_t m){
    mult_factor_vec.reserve(m);
    recursed_betas.reserve((m + 1) * m / 2);
    recursed_shifts.reserve(m);
    dimension = m;
  }
};

template <typename FieldT> 
struct PreComputedValues_Level2{
  std::vector<std::vector<FieldT>> sums_vec;
  std::vector<std::vector<FieldT>> betai_vec;
  size_t dimension;

  PreComputedValues_Level2(const size_t m){
    sums_vec.reserve(m);
    betai_vec.reserve(m);
    dimension = m;
  }
};

template <typename FieldT> 
struct PreComputedValues_CO{
  std::vector<FieldT> recursed_shifts;
  std::vector<FieldT> beta_sums_vec;
  size_t dimension;

  PreComputedValues_CO(const size_t m){
    recursed_shifts.reserve(m);
    dimension = m;
  }
};

template<typename FieldT>
PreComputedValues_Level1<FieldT> 
pre_computation_lvl1(const libiop::affine_subspace<FieldT> &domain);

template<typename FieldT>
PreComputedValues_Level2<FieldT> 
pre_computation_lvl2(const libiop::affine_subspace<FieldT> &domain);

template<typename FieldT>
PreComputedValues_CO<FieldT> 
pre_computation_co(const libiop::affine_subspace<FieldT> &domain);

template<typename FieldT>
std::vector<FieldT> 
additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                 const PreComputedValues_Level2<FieldT> &values);

template<typename FieldT>
std::vector<FieldT> 
additive_FFT(const std::vector<FieldT> &poly_coeffs,
                                 const PreComputedValues_Level1<FieldT> &values);

template<typename FieldT>
std::vector<FieldT> additive_FFT_CO(const std::vector<FieldT> &poly_coeffs,
                                 const PreComputedValues_CO<FieldT> &values);

template <typename FieldT>
std::vector<FieldT>
additive_FFT_CO(const std::vector<FieldT> &poly_coeffs,
                const libiop::affine_subspace<FieldT> &domain);

template <typename FieldT>
std::vector<FieldT>
additive_IFFT(const std::vector<FieldT> &evals,
              const libiop::affine_subspace<FieldT> &domain);

} // namespace gao

#include "Gao/fft.tcc"

#endif // ADDITIVE_FFT_GAO_HPP_
