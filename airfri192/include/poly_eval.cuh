#ifndef POLY_EVAL_HPP
#define POLY_EVAL_HPP

#include <vector>
#include <cstdint>
#include <libff/algebra/fields/binary/gf256.hpp>

// C++ interface for polynomial evaluation (wraps CUDA implementation)

namespace poly_eval {

/**
 * Evaluate polynomial at all points in domain (CPU-only version for now)
 * Later can be replaced with GPU-accelerated version
 * 
 * @param poly_coeffs: Polynomial coefficients (degree 0 to n-1)
 * @param domain_elements: Points to evaluate at
 * @return Codeword (evaluations at each domain point)
 */
std::vector<libff::gf256> evaluate_polynomial_on_domain(
    const std::vector<libff::gf256>& poly_coeffs,
    const std::vector<libff::gf256>& domain_elements
);

#ifdef __CUDACC__
// CUDA-accelerated version (only compiled when building with nvcc)
// This will be the GPU implementation
extern "C" void parallel_poly_eval_cuda(
    uint64_t **codeword, 
    uint64_t *poly_coeffs, 
    int poly_len, 
    uint64_t **domain_elements, 
    size_t field_words, 
    size_t initial_domain_length
);
#endif

} // namespace poly_eval

#endif // POLY_EVAL_HPP
