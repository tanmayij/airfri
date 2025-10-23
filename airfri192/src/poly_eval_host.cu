#include "poly_eval.hpp"
#include <iostream>

namespace poly_eval {

// CPU-only implementation of polynomial evaluation
// Evaluates poly at each domain element: result[i] = poly(domain[i])
std::vector<libff::gf256> evaluate_polynomial_on_domain(
    const std::vector<libff::gf256>& poly_coeffs,
    const std::vector<libff::gf256>& domain_elements
) {
    using FieldT = libff::gf256;
    
    size_t domain_size = domain_elements.size();
    size_t poly_len = poly_coeffs.size();
    std::vector<FieldT> codeword(domain_size);
    
    std::cout << "Evaluating polynomial (degree " << poly_len - 1 
              << ") on domain of size " << domain_size << std::endl;
    
    // For each domain element
    for (size_t i = 0; i < domain_size; ++i) {
        FieldT result = FieldT::zero();
        FieldT x_power = FieldT::one();
        const FieldT& x = domain_elements[i];
        
        // Horner's method: result = a0 + x*(a1 + x*(a2 + x*(...)))
        // Or direct: result = sum(a_j * x^j)
        for (size_t j = 0; j < poly_len; ++j) {
            result += poly_coeffs[j] * x_power;
            x_power *= x;
        }
        
        codeword[i] = result;
    }
    
    std::cout << "Polynomial evaluation complete." << std::endl;
    return codeword;
}

} // namespace poly_eval
