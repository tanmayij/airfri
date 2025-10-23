// Example integration of poly_eval into test_fri.cu

// At the top of test_fri.cu, add:
#include "poly_eval.hpp"

// Then you can use it in two ways:

// === Option 1: CPU polynomial evaluation (alternative to FFT) ===
// This directly evaluates poly(x) at each domain point
void example_cpu_poly_eval() {
    typedef libff::gf256 FieldT;
    
    // Your polynomial coefficients
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(32768); // degree 32767
    
    // Your domain (from Cantor basis)
    std::vector<FieldT> basis = cantor_basis<FieldT>(20);
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    
    // Evaluate polynomial on domain using CPU
    std::vector<FieldT> codeword = poly_eval::evaluate_polynomial_on_domain(
        poly_coeffs,
        domain.all_elements()
    );
    
    cout << "Codeword size: " << codeword.size() << endl;
}

// === Option 2: Keep using additive FFT (RECOMMENDED) ===
// Your current approach is already optimal!
void example_fft_approach() {
    typedef libff::gf256 FieldT;
    const size_t m = 20;
    
    // Generate random polynomial
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    
    // Use Cantor additive FFT (this is what you're already doing)
    int num_threads = std::thread::hardware_concurrency() > 0 ? 
                      std::thread::hardware_concurrency() - 1 : 31;
    std::vector<FieldT> codeword = cantor::additive_FFT(poly_coeffs, m, num_threads);
    
    cout << "Codeword size: " << codeword.size() << endl;
}

// === Note on GPU acceleration ===
// The .cu file (poly_eval.cu) is ready for GPU implementation, but:
// 1. You need to port field operations from FRI192/include/field.cuh
// 2. FFT is generally faster than direct evaluation for large domains
// 3. For now, stick with cantor::additive_FFT - it's already optimized

// Your current test_fri.cu is already doing the right thing!
// The poly_eval module is there if you need direct evaluation later.
