#include "libiop/algebra/subspace.hpp"
#include "libiop/algebra/utils.hpp"
#include "libiop/fft.hpp"
#include <algorithm>
#include <cstddef>
#include <libff/algebra/fields/binary/gf256.hpp>
#include <libff/algebra/fields/binary/gf128.hpp>
#include <libff/algebra/fields/binary/gf64.hpp>
#include <libff/algebra/fields/binary/gf32.hpp>
#include <libff/common/utils.hpp>
#include <ostream>
#include <vector>
#include <chrono>
#include <cmath>   // for std::sqrt

#include "utils/utils.hpp"
#include "Cantor/fft.hpp"
#include "Gao/fft.hpp"

template <typename T>
bool check_equal(const std::vector<T> &v1, const std::vector<T> &v2)
{
    for (size_t i = 0; i < v1.size(); ++i)
    {
        if (v1[i] != v2[i])
        {
            return false;
        }
    }
    return true;
}

/* This test uses different domains; hence, the timing report might be more accurate due to minimizing CPU's caching*/
void Cantor_FFT_Test(){
    typedef libff::gf256 FieldT;
    size_t m = 15;
    size_t N_test = 100;
    std::cout << "m = " << m << ", N_test = " << N_test << ", F = GF(2^"<<FieldT::extension_degree()<<")" << std::endl;
    std::vector<double> durations_Cantor(N_test);
    std::vector<double> durations_Gao(N_test);

    double Cantor_sum = 0;
    double Gao_sum = 0;

    std::vector<FieldT> basis(cantor_basis<FieldT>(m));

    for (size_t i=0; i<N_test; ++i){
        std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);

        // Cantor test
        libiop::field_subset<FieldT> cantor_domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
        auto start = std::chrono::high_resolution_clock::now();
        const std::vector<FieldT> cantor_result = cantor::additive_FFT<FieldT>(poly_coeffs, cantor_domain.subspace());
        auto stop = std::chrono::high_resolution_clock::now();
        durations_Cantor[i] = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        Cantor_sum += durations_Cantor[i];

        // Gao test
        // libiop::field_subset<FieldT> gao_domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
        libiop::field_subset<FieldT> gao_domain = cantor_domain;
        start = std::chrono::high_resolution_clock::now();
        const std::vector<FieldT> additive_result = libiop::additive_FFT<FieldT>(poly_coeffs, gao_domain.subspace());    // my_print_vector<FieldT>(cantor_result);
        stop = std::chrono::high_resolution_clock::now();
        durations_Gao[i] = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start).count();
        Gao_sum += durations_Gao[i];

        std::cout << "Equality check: " << (check_equal<FieldT>(additive_result, cantor_result) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m")  << std::endl;

    }
    
    double cantor_variance_sum = 0.0;
    for (double duration : durations_Cantor) {
        cantor_variance_sum += (duration - (Cantor_sum / N_test)) * (duration - (Cantor_sum / N_test));
    }

    double gao_variance_sum = 0.0;
    for (double duration : durations_Gao) {
        gao_variance_sum += (duration - (Gao_sum / N_test)) * (duration - (Gao_sum / N_test));
    }

    std::cout << "Cantor's average duration: " << Cantor_sum / N_test << " ms \t Standard deviations: " << std::sqrt(cantor_variance_sum / N_test) << std::endl;
    std::cout << "Gao's    average duration: " << Gao_sum    / N_test << " ms \t Standard deviations: " << std::sqrt(gao_variance_sum    / N_test) << std::endl;

}

/* This test is primarely for correctness check. We use the same domain; hence, the timing report might not be accurate due to CPU's caching*/
void Gao_CO_FFT_Test(){
    typedef libff::gf128 FieldT;
    size_t m = 20;
    std::cout << "m = " << m << ", Start testing!\n";

    std::vector<FieldT> basis(cantor_basis<FieldT>(m));
    std::reverse(basis.begin(), basis.end());

    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    gao::PreComputedValues_CO<FieldT> values = gao::pre_computation_co<FieldT>(domain.subspace());
    
    std::chrono::time_point<std::chrono::high_resolution_clock> start, stop; 

    std::cout << "Entering gao_co:\n";
    start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> gao_co_result = gao::additive_FFT_CO<FieldT>(poly_coeffs, domain.subspace());
    stop = std::chrono::high_resolution_clock::now();
    auto duration_gao_co = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "Entering gao_co precomputed:\n";
    start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> gao_co_result_precmp = gao::additive_FFT_CO<FieldT>(poly_coeffs, values);
    stop = std::chrono::high_resolution_clock::now();
    auto duration_gao_co_precmp = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);


    std::cout << "Entering libop:\n";
    start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> additive_result = libiop::additive_FFT<FieldT>(poly_coeffs, domain.subspace()); 
    stop = std::chrono::high_resolution_clock::now();
    auto duration_Gao = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "Equality check 1: " << (check_equal<FieldT>(additive_result, gao_co_result) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m") << std::endl;
    std::cout << "Equality check 2: " << (check_equal<FieldT>(additive_result, gao_co_result_precmp) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m") << std::endl;
    
    std::cout << "Gao_CO's Duration (pre-computed): " << duration_gao_co_precmp.count() << " ms" << std::endl;
    std::cout << "Gao_CO's Duration: " << duration_gao_co.count() << " ms" << std::endl;
    std::cout << "Gao's Duration: " << duration_Gao.count() << " ms" << std::endl;

}

/* This test is primarely for correctness check. We use the same domain; hence, the timing report might not be accurate due to CPU's caching*/
void Cantor_FFT_PreComputation_Test(){
    typedef libff::gf128 FieldT;
    size_t m = 15;
    std::cout << "m = " << m << ", Start testing!\n";

    size_t N_test = 10;

    for (size_t i=0; i<N_test; ++i){
    std::cout<<"test "<< i+1 << std::endl;
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    
    cantor::PreComputedValues<FieldT> values = cantor::pre_computation(domain.subspace());

    auto start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> cantor_precmp_result = cantor::additive_FFT<FieldT>(poly_coeffs, values);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration_Cantor_precmp = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> cantor_result = cantor::additive_FFT<FieldT>(poly_coeffs, domain.subspace());
    stop = std::chrono::high_resolution_clock::now();
    auto duration_Cantor = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    std::cout << "Equality check: " << (check_equal<FieldT>(cantor_precmp_result, cantor_result) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m")  << std::endl;
    std::cout << "Cantor's Duration (pre computed): " << duration_Cantor_precmp.count()<< " ms" << std::endl;
    std::cout << "Cantor's Duration: " << duration_Cantor.count()<< " ms" << std::endl;

    }
}

/* This test is primarely for correctness check. We use the same domain; hence, the timing report might not be accurate due to CPU's caching*/
void Gao_FFT_PreComputation_Test(){
    typedef libff::gf256 FieldT;
    size_t m = 18;
    std::cout << "m = " << m << ", Start testing!\n";

    libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);

    auto start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> gao_result = libiop::additive_FFT<FieldT>(poly_coeffs, domain.subspace());
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration_Gao = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    gao::PreComputedValues_Level1<FieldT> values_lvl1 = gao::pre_computation_lvl1<FieldT>(domain.subspace());
    start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> gao_precmp_result_lvl1 = gao::additive_FFT<FieldT>(poly_coeffs, values_lvl1);
    stop = std::chrono::high_resolution_clock::now();
    auto duration_gao_precmp_lvl1 = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);

    gao::PreComputedValues_Level2<FieldT> values_lvl2 = gao::pre_computation_lvl2<FieldT>(domain.subspace());
    start = std::chrono::high_resolution_clock::now();
    const std::vector<FieldT> gao_precmp_result_lvl2 = gao::additive_FFT<FieldT>(poly_coeffs, values_lvl2);
    stop = std::chrono::high_resolution_clock::now();
    auto duration_gao_precmp_lvl2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);


    std::cout << "Level1 Equality check: " << (check_equal<FieldT>(gao_precmp_result_lvl1, gao_result) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m")  << std::endl;
    std::cout << "Level2 Equality check: " << (check_equal<FieldT>(gao_precmp_result_lvl2, gao_result) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m")  << std::endl;
    std::cout << "Gao's Duration (pre computed level1): " << duration_gao_precmp_lvl1.count()<< " ms" << std::endl;
    std::cout << "Gao's Duration (pre computed level2): " << duration_gao_precmp_lvl2.count()<< " ms" << std::endl;
    std::cout << "Gao's Duration: " << duration_Gao.count()<< " ms" << std::endl;

}

void Valgrid_libiop_test(){
    typedef libff::gf256 FieldT;
    size_t m = 20;
    std::cout << "m = " << m << ", Start testing!\n";

    //Domain Creation
    libiop::field_subset<FieldT> domain = libiop::field_subset<FieldT>(libiop::affine_subspace<FieldT>::random_affine_subspace(m));
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);

    const std::vector<FieldT> gao_result = libiop::additive_FFT<FieldT>(poly_coeffs, domain.subspace());
}

void Valgrid_cantor_test(){
    typedef libff::gf256 FieldT;
    size_t m = 10;
    std::cout << "m = " << m << ", Start testing!\n";

    //Cantor Special Basis Domain Creation
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};
    
    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);

    const std::vector<FieldT> cantor_result = cantor::additive_FFT<FieldT>(poly_coeffs, domain.subspace());
}


void Valgrid_cantorPC_test(){
    typedef libff::gf256 FieldT;
    size_t m = 15;
    std::cout << "m = " << m << ", Start testing!\n";

    //Cantor Special Basis Domain Creation
    std::vector<FieldT> basis(cantor_basis<FieldT>(m));
    libiop::field_subset<FieldT> domain{libiop::affine_subspace<FieldT>(basis, FieldT::random_element())};

    std::vector<FieldT> poly_coeffs = libiop::random_vector<FieldT>(1ull << m);
    cantor::PreComputedValues<FieldT> values = cantor::pre_computation(domain.subspace());

    const std::vector<FieldT> cantor_result = cantor::additive_FFT<FieldT>(poly_coeffs, values);
    const std::vector<FieldT> poly_coeffs_computed = cantor::additive_IFFT<FieldT>(cantor_result, values);

    std::cout << "Equality check: " << (check_equal<FieldT>(poly_coeffs_computed, poly_coeffs) ? "\033[1;32mPass\033[0m" : "\033[1;31mFail\033[0m")  << std::endl;


}


void Pre_Compute_Cantor_basis(){
    typedef libff::gf192 FieldT;
    int N = 32;

    std::vector<FieldT> basis = cantor_basis<FieldT>(N);
    std::cout <<"cantor_in_gf"<<FieldT::extension_degree()<<"["<<N<<"]["<<FieldT::extension_degree()/64<<"] = {\n"; 
    for (int i = 0; i < N; i++){
        std::vector<uint64_t> basis_in_words = basis[i].to_words();
    std::cout <<"{ "; 
        for (int j = 0; j < basis_in_words.size() - 1; j++){
            std::cout<< basis_in_words[j] << ", ";
        }
        std::cout<< basis_in_words.back() << " }, \n";
    }
    std::cout <<"};\n"; 


}

int main()
{
    Pre_Compute_Cantor_basis();
    // Cantor_FFT_Test();
    // Gao_CO_FFT_Test();
    // Cantor_FFT_PreComputation_Test();
    // Gao_FFT_PreComputation_Test();

    // Valgrid_libiop_test();
    // Valgrid_cantorPC_test();

    return 0;
}
