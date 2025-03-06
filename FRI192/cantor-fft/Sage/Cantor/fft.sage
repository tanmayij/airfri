import math

def fft(g_coeffs, m, a, ext_degree, affine_shift=0):
    """
    Computes the Fast Fourier Transform (FFT) of the polynomial g(x) using Cantor's algorithm.

    Parameters:
    g_coeffs (list): A list of coefficients representing the polynomial g(x), ordered from the constant term 
                     to the highest degree term. After the FFT, this list will hold the evaluations of g(x) 
                     at points in the Cantor's basis.    
    m (int): An upper bound on the degree of g(x), ensuring that deg(g) < 2^m. This controls the size of the FFT.
    a (element of GF(2^ext_degree)): The primitive element of the finite field GF(2^ext_degree), which is used 
                                     to define the FFT points. It serves as the root of the irreducible polynomial 
                                     defining the finite field.
    ext_degree (int): The extension degree of GF(2), such that the finite field is GF(2^ext_degree). This controls 
                      the field in which the FFT computations are performed. In Cantor's algorithm, ext_degree 
                      must be a power of two.

    * The FFT is computed in place, so g_coeffs is modified to contain the evaluations of g(x) over Cantor's basis (W).
      Initially, g_coeffs represents the coefficients of the polynomial g(x), but by the end of the algorithm, it will 
      hold the values of g(x) evaluated at specific points from    
    """
    g_coeffs += [0]*(2**(m)-len(g_coeffs)) 
    _, nz_hdt_S, table = fft_precmp(a, m, ext_degree, affine_shift)
    fft_no_precmp_parallel(g_coeffs, m, nz_hdt_S, table)

def fft_no_precmp(g_coeffs, m, nz_hdt_S, table):
    """
    Performs the core computations of the Fast Fourier Transform (FFT) algorithm, excluding the precomputation of 
    S(x) polynomials and the shift table. This function assumes that these precomputations have already been done.

    Parameters:
    g_coeffs (list): A list of coefficients representing the polynomial g(x). This list will be modified in place 
                     during the FFT and will ultimately contain the evaluations of g(x) over the points in Cantor's basis.
    m (int): deg(g(x)) < 2^m. This corresponds to the number of stages in the FFT.
    nz_hdt_S (list[list]): A list of lists, where each inner list contains the indices of the non-zero coefficients 
                           in the polynomial s_r(x) for each round r. These indices are ordered from the highest 
                           degree term (hdt) to the lowest degree term.
    table (list[list]): A table of precomputed shifts used for tail module computations.

    * The result of the FFT is stored directly in g_coeffs, which is updated in place.
    """
    g_coeffs += [0]*(2**(m)-len(g_coeffs))
    # print(len(table), "here")
    input_size = len(g_coeffs)
    n_modules = 1
    # table = [item for sublist in table for item in sublist]
    cnt = 0
    for r in range(0,m):
        offset = 0
        for i in range(0, n_modules):
            tail_module(g_coeffs, nz_hdt_S[r], input_size, offset, table[cnt])
            cnt += 1
            offset += input_size
        input_size >>= 1
        n_modules  <<= 1

def fft_no_precmp_parallel(g_coeffs, m, nz_hdt_S, table):
    g_coeffs += [0]*(2**(m)-len(g_coeffs))
    input_size = len(g_coeffs)
    n_modules = 1
    cnt = 0
    for r in range(0,m):
        offset = 0
        half_input_size = input_size >> 1
        for i in range(0, n_modules):
            mult_factor = table[cnt]
            cnt += 1
            offset2 = offset + half_input_size
            for k in range(offset2+half_input_size -1, offset2 -1, -1):
                g_k = g_coeffs[k]
                for nz in nz_hdt_S[r][:-1]: 
                    g_coeffs[k - nz] += g_k
                g_coeffs[k-half_input_size] += g_k * mult_factor
            for j in range(half_input_size):
                g_coeffs[offset2+j] += g_coeffs[offset+j]  
            offset += input_size
        input_size = half_input_size
        n_modules  <<= 1




def fft_precmp(a, m, ext_degree, affine_shift=0):
    """
    Performs the precomputations required for Cantor's FFT algorithm. These precomputations include
    (1) constructing Cantor's basis, (2) computing the polynomials s_r(x) for each FFT round and determining 
    the corresponding indices of their non-zero coefficients, and (3) generating a table of precomputed shifts 
    for use in the tail module computations.

    Parameters: 
    a (element of GF(2^ext_degree)): The primitive element of the finite field GF(2^ext_degree). It is the root of the 
                                     irreducible polynomial that defines the field.
    m (int): deg(g(x)) < 2^m. This determines the size of the FFT and the number of rounds.
    ext_degree (int): The extension degree of the field GF(2), such that the field is GF(2^ext_degree). In Cantor's algorithm, 
                      ext_degree must be a power of two to compute Canto's basis.

    Returns:
    W (list): Cantor's basis, which defines the set of points where the polynomial will be evaluated. 
    S (list[list]): A list of polynomials s_r(x) for each FFT round r. Each polynomial s_r(x) is represented as a list of 
                    its coefficients. These polynomials are used in the modular reductions at each stage of the FFT.
    nz_hdt_S (list[list]): A list of lists where each inner list contains the indices of the non-zero coefficients 
                           in the polynomial s_r(x) for each round r, ordered from the highest degree term (hdt) to the lowest.
    table (list[list]): A precomputed table of shifts used in tail module computations during each FFT round.

    """
    W = fast_initial_basis_computation(a, m, ext_degree)
    nz_hdt_S, affine_shift_list = S_function_computation_2(m, affine_shift)
    table = fast_S_shifts_table_generator(W, affine_shift_list)
    # print(table)
    # print(affine_shift_list)
    # quit()
    return W, nz_hdt_S, table



def S_function_computation(m):
    """
    Computes the s_r(x) polynomials for each round r where 0 <= r < m.

    Parameters: 
    m (int): Determines the number of rounds (r), where r ranges from 0 to m-1.

    Returns: 
    nz_hdt_S (list[list]): A list of lists where each inner list contains the indices of the non-zero coefficients 
                           in the polynomial s_r(x) for each round r, ordered from the highest degree term (hdt) to the lowest.
    """
    # nz_hdt_S: Non zero coeffs indices in S ordered from the highest degree term (hdt) to the lowest.
    # Hence, for example  0 corresponds to the highest degree term in s_r(x).
    nz_hdt_S = [[]] * (m)                   
    for r in range(m):
        # S[r] = [0] * (2**r + 1) # deg(S_r) = 2^r because 'C(r,r) = 1' 
        nz_hdt_S[r] = []
        for i in range(r+1): # S_r(t) = \sum_{i=0}^{r} C(r,i)*(t^{2^i}). Therefore, the loop should include i = r
            # S[r][2**i] = math.comb(r,i) % 2 # Coefficients are ordered from the constant term up to the highest degree term.
            if (math.comb(r,i) % 2): nz_hdt_S[r].append((1<<r) - (1<<i)) # Non-zero coefficients indices in S ordered from highest degree term (hdt)       
    nz_hdt_S.reverse()  # Reverse the list because the FFT algorithm processes the polynomials starting from S_{m-1} 
                        # down to S_0.
    # S.reverse()
    return nz_hdt_S

def S_function_computation_2(m, affine_shift=0):
    """
    Computes the s_r(x) polynomials for each round r where 0 <= r < m. This implementation avoids the use of math.comb 
    and instead utilizes bitwise operations to determine the parity of binomial coefficients C(r, i).
    
    In Sage, this implementation does not show improved performance compared to the original `S_function_computation` 
    because math.comb is implemented in C, which makes it highly efficient. However, in C++ we expect this implementation 
    to offer better performance due to the efficiency of bitwise operations.

    Parameters: 
    m (int): Determines the number of rounds (r), where r ranges from 0 to m-1.

    Returns: 
    nz_hdt_S (list[list]): A list of lists where each inner list contains the indices of the non-zero coefficients 
                           in the polynomial s_r(x) for each round r, ordered from the highest degree term (hdt) to the lowest.
    """
    # nz_hdt_S: Non zero coeffs indices in S. Note that index 0 is the coefficient of the highest degree term in s_r(x),
    affine_shift_list = [affine_shift] * m
    nz_hdt_S = [[]] * (m)                   
    for r in range(m-1):
        S_index = m + ~r
        nz_hdt_S[r] = [(1<<S_index)-1]
        for i in range(1,S_index): # S_r(t) = \sum_{i=0}^{r} C(r,i)*(t^{2^i}). 
                             # C(r,0) = 1, so we assigned that before loop.
                             # Also, C(r,r) = 1, so we will append that after loop.
            is_odd = (S_index & True) or (not(i & True))
            ii = i >> 1
            rr = S_index >> 1
            while(is_odd and ii>0):
                is_odd = (rr & True) or (not(ii & True))
                rr >>= 1
                ii >>= 1
            if (is_odd): 
                nz_hdt_S[r].append((1<<S_index) - (1<<i)) # Non-zero coefficients indices in S ordered from highest degree term (hdt)       
                affine_shift_list[r] += affine_shift**(1<<i)
        nz_hdt_S[r].append(0)
        affine_shift_list[r] += affine_shift**(1<<S_index)
    return nz_hdt_S, affine_shift_list

def S_function_computation_3(m):
    """
    Computes the s_r(x) polynomials for each round r where 0 <= r < m using the recursive relation:
    S_{r+1}(x) = S(S_r(x)) = S_r^2(x) + S_r(x).

    Parameters: 
    m (int): Determines the number of rounds (r), where r ranges from 0 to m-1.

    Returns: 
    nz_hdt_S (list[list]): A list of lists where each inner list contains the indices of the non-zero coefficients 
                           in the polynomial s_r(x) for each round r, ordered from the highest degree term (hdt) to the lowest.
    """
    S = [[]] * (m)
    nz_S = [[]] * (m) 
    nz_hdt_S = [[]] * (m) 
    S[0] = [0,1] # S_0 = x
    nz_S[0] = [1]
    nz_hdt_S[0] = [0]
    for r in range(1,m):
        S[r] = ([0] * ((1<<(r)) + 1)) # deg(S_r) = 2^r hence it has at most 2^r + 1 coefficients
        nz_hdt_S[r] = []
        nz_S[r] = []
        for i in nz_S[r-1]:
            S[r][2*i] = S[r-1][i] 
            S[r][i] ^^= S[r-1][i] 

        for i in range((1<<r)+1):
            if (S[r][i] == 1): 
                nz_hdt_S[r].append((1<<r) - i)
                nz_S[r].append(i)

    nz_hdt_S.reverse()
    S.reverse()
    return nz_hdt_S
    
        
def S_shifts_table_generator(S, W):
    """
    Computes the shift values required for the tail_module function. In each branch of each round, 
    the tail_module function requires the value of S_r(y_p + y_q), where y_p and y_q are either one of the basis elements in W
    or are zero. This function pre-computes all the non-zero values by evaluating them in the S_r polynomial.
 
    Parameters: 
    S (list[list]): A list of lists, where each sublist contains the coefficients of the polynomial S_r for round r.
    W (list): A list of basis elements.

    Returns: 
    table (list[list]): A list of lists, where each inner list contains the precomputed shift values for the tail_module function 
                        for each branch in each round.
    """
    table = []
    for i in range(len(W)):
        row = []
        for j in range((2**(i+1))):
            if i == 0:
                shift = 0
            else:
                shift = table[-1][j>>1]
            if j % 2:
                shift += W[~i]          # ~i points to the i-th last element in the list:  ~ is the bitwise NOT 
                                        #where ~i = -i-1. bitwise NOT is used because its time complexity is O(1)
            row += [shift]
        table += [row]    
    for i in range(len(table)):
        for j in range(len(table[i])):
            table[i][j] = eval_s_at_x(S[i], table[i][j])
        table[i] = table[i][::2]
        table[i] = table[i][1:]

    return table

def fast_S_shifts_table_generator(W, affine_shift_list = None):
    """
    Computes the shift values required for the tail_module function. In each branch of each round, 
    the tail_module function requires the value of S_r(y_p + y_q), where y_p and y_q are either one of the basis elements in W
    or are zero. This function optimizes the computation by not explicitly evaluating y_p and y_q in the S_r polynomials.
    Instead, it uses the fact that S_r(y_{r+k}) = y_k, where y_k in W.
 
    Parameters: 
    S (list[list]): A list of lists, where each sublist contains the coefficients of the polynomial S_r for round r.
    W (list): A list of basis elements.

    Returns: 
    table (list[list]): A list of lists, where each inner list contains the precomputed shift values for the tail_module function 
                        for each branch in each round.
    """
    m = len(W)
    table = [0] * ((1<<m)-1)
    if affine_shift_list == None:
        affine_shift_list = [0] * len(W)
    cnt = 0
    for r in range(0, len(W)): # r is the row number. The first row is skipped as it only has a head module.
        # row = [affine_shift_list[r]] * ((1<<(r)) )  # 1<<(r) = 2^r which is the number of modlues in the row r. 
                                    # The first module is skipped as it is a head module.     
        table[cnt] = affine_shift_list[r]
        cnt += 1
        for module in range(1, 1<<(r)): 
            table[cnt] += affine_shift_list[r]
            for i in range(r):
                if module & (1<<i): # if the i-th bit in the r-bit module number is set, add W[i+1] to the shift value.
                    table[cnt] += W[i+1]
            cnt +=1
    return table


def eval_s_at_x(s, x):
    """
    Evaluates the polynomial s at a given point x in GF(2^ext_degree), where the coefficients of s are in GF(2).
    
    Parameters:
    s (list):  A list of coefficients for the polynomial s, where each coefficient is in GF(2). 
               The i-th element of the list corresponds to the coefficient of x^i.
    x (element of GF(2^ext_degree)): The point at which the polynomial is evaluated. 
                                     It is an element of the finite field GF(2^ext_degree).

    Returns:
    result (element of GF(2^ext_degree)): The result of the evaluation, i.e., s(x). 
                                          This is the value of the polynomial at the point x.
    """
    result = 0
    for i in range(len(s)):
        result += s[i] * (x ** i)
    return result

def divide(coeffs, nz_hdt_S, input_size, offset, multـfactor=1):
    """
    Computes the quotient q(x) = g(x) / s(x), where the dividend polynomial g(x) is represented 
    by the coefficients in coeffs[offset:offset+input_size], and the divisor polynomial s(x) 
    is represented by the non-zero indices in nz_hdt_S. s(x) is in GF(2).

    Parameters:
    coeffs (list): List of coefficients representing the polynomials. This is modified in place 
                   to store the remainder after division.
    nz_hdt_S (list): List of indices representing the non-zero terms of the divisor polynomial s(x), 
                     where the coefficients of s(x) are in GF(2).
    input_size (int): Number of coefficients in the dividend polynomial g(x).
    offset (int): Starting index in 'coeffs' where the dividend polynomial g(x) begins.

    Returns:
    q (list): Coefficients of the quotient polynomial q(x).
    """
    # Initialize the quotient polynomial q(x), whose degree is at most half the degree of the dividend g(x).
    # q = [0] * (input_size>>1)       # The degree of q(x) is derived as deg(g(x)) - deg(s(x)).
                                    # For input size 2^(r+1), we have deg(q) <= 2^r - 1 = input_size // 2 - 1.        
    half_input_size = input_size>>1
    offset += half_input_size

    # Process each term in the quotient by iterating over the coefficients in reverse order.
    for i in reversed(range(offset,offset+half_input_size)):   # The reasorown to use 'reversed()' is that coefficients are ordered from 
                                                                        # the constant term up to the highest degree term. 
                                                                        # Therefore, the last element is the highest degree coefficient.
        # q[~q_ind] = coeffs[i]       # Using ~q_ind accesses elements from the end.

        for nz in nz_hdt_S[:-1]:    # Exclude the highest degree term of s(x) (the last element of nz_hdt_S),
                                    # because it cancels out with coeffs[i].
            coeffs[i - nz] += coeffs[i]
        coeffs[i-half_input_size] += coeffs[i] * multـfactor
        # coeffs[i] += coeffs[i-len(q)]

        # q_ind += 1

def head_module(coeffs, nz_hdt_S, input_size):
    """
    The head_module is the first module in each row. It computes on the first input_size coefficients of coeffs.

    Parameters:
    coeffs (list):  A list of coefficients representing the polynomial g(x) at the current round
    nz_hdt_S (list): A list of indices representing the non-zero terms of the polynomial s_r(x), ordered from the highest degree term to the lowest. 
                     This list is used as input to the modules (both head or tail) in round r
    input_size (int): The number of coefficients in coeffs that are provided as input to this module.

    * This module modifies coeffs in place.
    """
    q = divide(coeffs, nz_hdt_S, input_size, offset=0)
    for i in range(len(q)):
        coeffs[i+len(q)] = coeffs[i] + q[i]

def tail_module(coeffs, nz_hdt_S, input_size, offset, s_shift1):
    """
    The tail_module is used as modules after the head_module in each row. It operations on a subset of coefficients, starting from a given offset.

    Parameters:
    coeffs (list):  A list of coefficients representing the polynomial g(x) at the current round
    nz_hdt_S (list): A list of indices representing the non-zero terms of the polynomial s_r(x), ordered from the highest degree term to the lowest. 
                     This list is used as input to the head module in round r
    input_size (int): The number of coefficients in coeffs that are provided as input to this module.
    offset (int): The starting index from which this module operates within coeffs
    s_shift1 (int): A shift factor applied to the result of the division

    * This module modifies coeffs in place.
    """
    divide(coeffs, nz_hdt_S, input_size, offset, s_shift1)
    # if (q != coeffs[offset+len(q):offset+2*len(q)]):
    #     print ("q:", q)
    #     print (f"coeffs[{offset}+{len(q)}:{offset}+2*{len(q)}:", coeffs[offset+len(q):offset+2*len(q)])
    #     quit()
    half_input_size = input_size >> 1
    for i in range(half_input_size):
        # coeffs[offset+i] += coeffs[offset+i+len(q)] * s_shift1
        coeffs[offset+i+half_input_size] += coeffs[offset+i]  
        pass