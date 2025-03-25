
# FRI192

## Overview
FRI192 is an implementation of the Fast Reed-Solomon Interactive Oracle Proof of Proximity (FRI) protocol designed to work with CUDA. The FRI protocol is used in cryptographic proof systems for low-degree testing of polynomials over finite fields. This project leverages GPU acceleration to improve performance and handle large inputs efficiently.

---

## Requirements
### Hardware
- **NVIDIA GPU** with CUDA support (CUDA Compute Capability 6.0+)

### Software
- **CUDA Toolkit** (tested with CUDA 11+)
- **OpenSSL** (for cryptographic operations)
- **GCC** (for compiling host code)
- **Make** (for build automation)

### Environment Setup
If you're using `conda`, you can set up the environment as follows:

```bash
conda create -n fri192_env python=3.9
conda activate fri192_env
conda install -c conda-forge cudatoolkit
conda install -c conda-forge openssl
conda install nvcc
```

---

## Build Instructions
### 1. Clone the repository
```bash
git clone <repository-url>
cd FRI192
```

### 2. Build the project
Use the `Makefile` provided to compile and link the project:
```bash
make
```

### 3. Run the program
To run the compiled binary:
```bash
make run
```

---

## File Structure
```
FRI192/
├── include/                 # Header files
├── src/                     # Source files
├── Makefile                 # Build automation
├── README.md                # Project documentation
```

---

## Code Overview
### Key Components:
| File                  | Description                                           |
|-----------------------|-------------------------------------------------------|
| **fri.cu**             | Main implementation of the FRI protocol              |
| **merkle.cu**          | Merkle tree construction and verification             |
| **commit-launch- merkle.cu**          | has GPU kernel functions               |
| **hash.cu**            | Cryptographic hashing utilities using SHA3            |
| **fft.cu**             | Fast Fourier Transform implementation                |
| **field.cu**           | Finite field arithmetic                              |
| **poly.cu**            | Polynomial interpolation and evaluation               |
| **fiat-shamir.cu**     | Fiat-Shamir transform for converting interactive proofs to non-interactive ones |
| **params.cu**          | Parameters and constants for the protocol             |

---

Not implemented yet
To debug CUDA memory issues:
```bash
cuda-gdb ./fri
```

Check GPU usage:
```bash
nvidia-smi
```

---

## Example Run
```bash
./fri
```
Example Output:
```
Colinearity check passed
Authentication Path:
Layer 0: aee28f0ca6e483fd 16785a4fdd12fc0d 80f790fbc4ceb528 b892ec4efd109db4 
...
Computed Merkle Root: 2f9731ece103d2be c939d3db9c450f95 17c12b335abffa62 5376d51a4a464fea 
Verification successful!
```

---

## Performance Notes
- Ensure that the GPU is set to **"Maximum Performance Mode"** using:
```bash
nvidia-smi -pm 1
```
- Disable GPU auto-throttling:
```bash
nvidia-smi -ac 5001,1410
```

---

## Benchmark Results
| Prover Time | Verifier Time | Degree     |
|-------------|---------------|------------|
|  0.670      | 0.166         | 0.01 sec   |


---

## Notes
Please contact me for the documentation of this implementation

---

## Contact
For questions or contributions, please open an issue or submit a pull request.

---
