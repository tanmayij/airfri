# Additive-FFT
Additive Fast Fourier Transforms over Finite Fields

## Overview

This repository contains both SageMath and C++ implementations related to Cantor and Gao-Mateer algorithms. The C++ implementation relies on the `libff` library, which is included as a Git submodule, so there's no need for a separate installation of `libff`. For the SageMath implementation, only Sage is required.



<!-- 
## To-do:
1. Implement Gao's general algorithm using Cantor's basis.
1. Implement Gao's special case algorithm
1. Implement the IFFT algorithm for Cantor's algorithm in Sage
1. Implement the IFFT algorithm for Gao's algorithm in Sage 
1. Design Gao's and Cantor's algorithms for the affine space
1. Implement Gao's and Cantor's algorithms for the affine space. I should change the current implement to general affine space implementations.

### Done:
 1. Make `S_shifts_table_generator` function faster in the pre computations in the Cantor's algorithm using the fact that $S_{i+k}(y_k) = y_i$.
 1. Write comments on the code
 1. Gao's FFT excluding the pre-computations


## Timings

### 1. $m = 14, GF(2^{256})$:
#### Cantor's FFT (average over 100 iterations)
- excludes pre-computation: 0.5657543277740479 s
- includes pre-computation: 0.7454671454429627 s
#### Gao's FFT (average over 10 iterations)
- includes pre-computation: 5.390587258338928 s


### 2. $m=11, GF(2^{256})$:
### Cantor's FFT (average over 100 iterations)
- Average Cantor's FFT time (excludes pre-computation): 0.06200750589370727 s
- Average Cantor's FFT time (includes pre-computation): 0.08558360099792481 s

### Gao's FFT over Cantor Basis (average over 100 iterations)
- Average Gao's FFT time (excludes pre-computation lvl2):                   0.11530064821243285 s
- Average Gao's FFT time (excludes pre-computation lvl1):                   0.1717957592010498 s
- Average Gao's FFT time (Full: includes pre-computation):                  0.17836191177368163 s
- Average Gao's Cantor optimized FFT time (excludes pre-computation lvl2):  0.10389098882675171 s
- Average Gao's Cantor optimized FFT time (Full: includes pre-computation): 0.15841137647628784 s


### 3. $m=11, GF(2^{256})$ Affine space: 
- Average direct evaluation time: 2.6934192204475402 s
- Average Cantor's FFT time (excludes pre-computation): 0.03471111297607422 s
- Average Cantor's FFT time (excludes pre-computation - parallel mode): 0.03284716844558716 s
- Average Cantor's FFT time (includes pre-computation - parallel mode): 0.04717706918716431 s
 -->
