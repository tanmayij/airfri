def span_basis(B, affine_shift=0):
    subset = [0] * 2**len(B)
    for i in range(len(subset)):
        subset[i] = an_element_in_basis(B, i) + affine_shift
    return subset

def polynomial_to_string(coeffs):
    string = f"{coeffs[0]}"
    for i, c in enumerate(coeffs[1:]):
        if c:
            if c==1:
                string += f" + x^{i+1}"
            else:
                string += f" + ({c})x^{i+1}"
    return string

def evaluate_polynomial(coeffs, eval_set):
    """
    Evaluates the polynomial represented by coeffs at each point in the evaluation set eval_set.

    Parameters:
    coeffs (list): The list of coefficients of the polynomial, ordered from the constant term up to the highest degree term.
    eval_set (list): The list of points at which to evaluate the polynomial.

    Returns:
    evaluations (list): A list of the polynomial evaluations at each point in eval_set.
    """
    evaluations = [0] * len(eval_set)
    for i in range(len(eval_set)):
        if eval_set[i] == 0:
            evaluations[i] = coeffs[0]
        else:
            x = 1
            for c in coeffs:
                evaluations[i] += c * x
                x *= eval_set[i]
    return evaluations


def an_element_in_basis(B, index):
    """
    Returns the element at the specified index in the basis set B.

    Parameters:
    B (list): The basis set.
    index (int): The index of the element to retrieve from the basis set B.

    Returns:
    result: The element at the specified index in B, s.t., B[index]
    """
    index_bin = Integer(index).binary()
    result = 0
    for i, bit in enumerate(reversed(index_bin)):
        result += int(bit) * B[i]
    return result


def initial_basis_computation(a, m):
    """
    Computes Cantor's basis non-optimally. The algorithm initializes W[0] = 1 and derives each subsequent W[i], such that 
    W[i]^2 + W[i] = W[i-1] for i >= 1. It performs a search in the field GF(2^ext_degree) to find each W[i]. 
    This is a brute-force approach and not an efficient method.
    
    Parameters: 
    a (element of GF(2^ext_degree)): The primitive element of the finite field GF(2^ext_degree). It is the root of the 
                                     irreducible polynomial that defines the field.
    m (int): deg(g(x)) < 2^m. This determines the size of the FFT and the number of rounds.
    
    Returns: 
    W (list): Cantor's basis, which defines the set of points where the polynomial will be evaluated. 
    """
    W = [1] * m
    for i in range(1, m):
        b = a
        while b**2 + b != W[i-1]:
            b *= a
        W[i] = b
    return W

def fast_initial_basis_computation(a, m, ext_degree):
    """
    Computes Cantor's basis optimally using the trace function over GF(2^ext_degree). The algorithm first searches for W[m-1] 
    such that tr(W[m-1]) = 1, where tr() is the trace function. Once W[m-1] is found, the algorithm recursively determines 
    W[i-1] for each i, such that W[i]^2 + W[i] = W[i-1].
     
    Parameters: 
    a (element of GF(2^ext_degree)): The primitive element of the finite field GF(2^ext_degree). It is the root of the 
                                     irreducible polynomial that defines the field.
    m (int): deg(g(x)) < 2^m. This determines the size of the FFT and the number of rounds.
    ext_degree (int): The extension degree of the field GF(2), such that the field is GF(2^ext_degree). In Cantor's algorithm, 
                      ext_degree must be a power of two to compute Canto's basis.

    Returns: 
    W (list): Cantor's basis, which defines the set of points where the polynomial will be evaluated. 
    """
    W = [1] * ext_degree
    b_m = a
    while(b_m.trace() != 1):
        b_m *= a
    W[-1] = b_m

    for i in reversed(range(1, ext_degree-1)):
        W[i] = W[i+1]**2 + W[i+1]
    return W[:m]