from time import time
import copy
import matplotlib.pyplot as plt



load('fft.sage')
load('../utils/utils.sage')

def test_fft(a, FF, ext_degree):
    N_tests = 1
    DIRECT_EVAUATION_TEST = True
    
    direct_eval_time = [0] * N_tests
    cantors_fft_no_precmp_time = [0] * N_tests
    cantors_fft_no_precmp_time_parallel = [0] * N_tests
    cantors_fft_with_precmp_time = [0] * N_tests

    m = 4
    print("Entering Pre-computation")
    affine_shift =  FF.random_element()
    # print("affine_shift:", affine_shift.to_integer())
    W, nz_hdt_S, table = fft_precmp(a, m, ext_degree, affine_shift)
    print("Finished Pre-computation")

    if DIRECT_EVAUATION_TEST:
        evaluation_set = span_basis(W, affine_shift)

    
    for iter in range(N_tests):
        g_coeffs = [FF.random_element() for i in range(2**m)]
        g_coeffs_copy = copy.deepcopy(g_coeffs)
        g_coeffs_copy_2 = copy.deepcopy(g_coeffs)


        if DIRECT_EVAUATION_TEST:
            print(iter, "Entering Direct Evaluation")
            start = time()
            evaluated_polynomial = evaluate_polynomial(g_coeffs, evaluation_set)
            direct_eval_time[iter] = time() - start

        print(iter, "Entering FFT excluding pre-computation")
        start = time()
        fft_no_precmp(g_coeffs, m, nz_hdt_S, table)        
        cantors_fft_no_precmp_time[iter] = time() - start

        if(DIRECT_EVAUATION_TEST and g_coeffs != evaluated_polynomial):
            print("Error: test failed for \"exclude pre-computation\"")
            exit()

        print(iter, "Entering FFT excluding pre-computation (parallel mode)")
        start = time()
        fft_no_precmp_parallel(g_coeffs_copy_2, m, nz_hdt_S, table)        
        cantors_fft_no_precmp_time_parallel[iter] = time() - start

        if(DIRECT_EVAUATION_TEST and g_coeffs_copy_2 != evaluated_polynomial):
            print("Error: test failed for \"exclude pre-computation (parallel mode)\"")
            exit()
        
        print(iter, "Entering FFT including pre-computation")
        start = time()
        fft(g_coeffs_copy, m, a, ext_degree, affine_shift)
        cantors_fft_with_precmp_time[iter] = time() - start

        if(g_coeffs_copy != g_coeffs):
            print("Error: test failed for \"includes pre-computation\"")
            exit()

    print("All tests passed")
    if DIRECT_EVAUATION_TEST:
        print("Average direct evaluation time:", sum(direct_eval_time)/N_tests, 's')
    print("Average Cantor's FFT time (excludes pre-computation):", sum(cantors_fft_no_precmp_time)/N_tests, 's')
    print("Average Cantor's FFT time (excludes pre-computation - parallel mode):", sum(cantors_fft_no_precmp_time_parallel)/N_tests, 's')
    print("Average Cantor's FFT time (includes pre-computation - parallel mode):", sum(cantors_fft_with_precmp_time)/N_tests, 's')

    test_numbers = list(range(N_tests))
    plt.figure(figsize=(12, 6))
    if DIRECT_EVAUATION_TEST:
        plt.plot(test_numbers, cantors_fft_with_precmp_time, marker='x', label='Direct')
        pass
    plt.plot(test_numbers, cantors_fft_with_precmp_time, marker='o', label='FFT Full - parallel')
    plt.plot(test_numbers, cantors_fft_no_precmp_time_parallel, marker='^', label='FFT precomputed - parallel')
    plt.plot(test_numbers, cantors_fft_no_precmp_time, marker='s', label='FFT precomputed')

    plt.xlabel('Test Number')
    plt.ylabel('Time (seconds)')  # Adjust the unit if necessary
    plt.title('FFT Timing Results Across Tests')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    return (sum(direct_eval_time)/N_tests, sum(cantors_fft_no_precmp_time)/N_tests, sum(cantors_fft_with_precmp_time)/N_tests)

    

def generate_a_map(a, dim):
    for i in range(2**dim):
        print(f"a^{i} = {a**i}")

if __name__ == "__main__":
    # F.<xx> = QQ[]
    # FF.<a> = GF(2**4, modulus= xx**4 + xx + 1)

    F.<x> = GF(2)[]
    ext_degree = 32
    irreducible_poly = F.irreducible_element(ext_degree)
    # irreducible_poly = x^32 + x^22 + x^2 + x^1 + 1 (to mach libiop for gf32)
    FF.<a> = GF(2**ext_degree, modulus=irreducible_poly)
    generate_a_map(a, 4)
    test_fft(a, FF, ext_degree)
    # fast_initial_basis_computation(a, 13)
