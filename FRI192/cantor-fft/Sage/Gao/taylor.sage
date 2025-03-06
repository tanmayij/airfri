
load('../utils/utils.sage')
import copy

def g_0_g_1_extraction(coeffs, input_size, b):
    coeffs_copy = [0] * len(coeffs)
    offset = 0
    for _ in range(b):
        coeffs_copy[offset:offset+input_size] = coeffs[offset:offset+input_size:2] + coeffs[offset+1:offset+input_size:2]
        coeffs [offset:offset+input_size] = coeffs_copy[offset:offset+input_size]
        offset += input_size
    # return coeffs_copy

def taylor_expansion_no_post_cmp(coeffs, input_size, input_offset):
    """
    Since, in the last part of the Taylor expansion, we need to arranges the even-indexed elements first, 
    followed by the odd-indexed elements. Hence, we have to copy
    """
    r = 0
    while input_size >> 2:
        offset = input_offset
        for _ in range(1<<r): 
            taylor_module(coeffs, input_size, offset) 
            offset += input_size
        r += 1
        input_size >>= 1

def taylor_expansion(coeffs, input_len):
    """
    Since, in the last part of the Taylor expansion, we need to arranges the even-indexed elements first, 
    followed by the odd-indexed elements. Hence, we have to copy
    """
    input_size = input_len
    r = 0
    while input_size >> 2:
        offset = 0
        for _ in range(1<<r): 
            taylor_module(coeffs, input_size, offset) 
            offset += input_size
        r += 1
        input_size >>= 1
    # print("g_0(x) = ", polynomial_to_string(coeffs[0::2]))
    # print("g_1(x) = ", polynomial_to_string(coeffs[1::2]))
    coeffs[:] = coeffs[0::2] + coeffs[1::2] # Concats g_0 and g_1
    return (coeffs) 
    

def taylor_module(coeffs, input_size, offset):
    chunk_size = input_size >> 2     # Divides the vector of coeffs into four chunks
    offset1 = offset  + chunk_size
    offset2 = offset1 + chunk_size
    offset3 = offset2 + chunk_size
    for i in range(chunk_size):
        coeffs[ offset2 + i ] += coeffs[ offset3 + i ]
        coeffs[ offset1 + i ] += coeffs[ offset2 + i ]
