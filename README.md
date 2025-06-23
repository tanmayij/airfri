# Test Plan and Commands

This document presents the test plan and benchmarking results for different implementations and optimisations of the FRI (Fast Reed–Solomon Interactive Oracle Proof of Proximity) protocol. It includes performance benchmarks under various settings and implementation phases across CPU and GPU variants, followed by essential compilation and execution commands.

---

## Test Plan for Various Optimisations (L3 Security Level)

| **S.No.** | **Category** | **Time (Seconds)** | **Finite Field** | **Degree** | **Expansion Factor** | **LDT Small Poly Degree** | **Initial Codeword Length** | **# Collinearity Tests** | **# Rounds** | **Merkle Tree Height** | **Hash Input → Output** | **Max Degree (Hardcoded)** |
|----------|--------------|--------------------|------------------|------------|-----------------------|----------------------------|------------------------------|--------------------------|-------------------------|------------------------|----------------------------|-----------------------------|
| 1 | FRI in PREON’s code | Avg: 143.93 (Prover)<br>Avg: 1.58 (Verifier)<br>52 runs | 2^256 | 4096 | 32 | 5 | 2^17 (131072) | 13 | 12 | 17 | 1024 → 256 | 15 |
| 2 | Thesis – C (no opt.) | Avg: 22.52 | 2^256 | 4096 | 32 | 5 | 2^17 (131072) | 13 | 12 | 17 | 1024 → 256 | 15 |
| 3 | Thesis – C (with precomp. inverses & Fiat-Shamir) | Avg: 17.95 | 2^256 | 4096 | 32 | 5 | 2^17 (131072) | 13 | 12 | 17 | 1024 → 256 | 15 |
| 4 | Thesis – C (with all optimisations incl. Merkle) | Avg: 13.47 | 2^256 | 4096 | 32 | 5 | 2^17 (131072) | 13 | 12 | 17 | 1024 → 256 | 15 |
| 5 | Thesis – GPU (precomp., parallel line comp.) | Avg: 9.36 | 2^256 | 4096 | 32 | 5 | 2^17 (131072) | 14 | 12 | 17 | 1024 → 256 | 15 |
| 6 | Thesis – GPU (with FRI Merkle Tree) | Avg: 1.49 | 2^256 | 4096 | 32 | 5 | 2^17 (131072) | 12 | 12 | 17 | 1024 → 256 | 15 |

---

## Test Plan for L5 Security Level

| **S.No.** | **Category** | **Finite Field** | **Degree** | **Expansion Factor** | **LDT Small Poly Degree** | **Max Degree (Hardcoded)** | **Initial Length** | **Final Length** | **# Rounds** | **Merkle Tree Height** | **Hash Input → Output** | **# Collinearity Tests** | **Time (Seconds)** |
|----------|--------------|------------------|------------|-----------------------|----------------------------|-----------------------------|--------------------|------------------|--------------|------------------------|----------------------------|--------------------------|--------------------|
| 1 | FRI in PREON’s code | 2^320 | 4096 | 32 | 5 | 15 | 2^17 | 2^5 | 12 | 17 | 1024 → 256 | 13 | Avg: 143.93 (Prover)<br>Avg: 1.58 (Verifier)<br>52 runs |
| 2 | Thesis – C (no opt.) | 2^320 | 4096 | 32 | 5 | 15 | 2^17 | 2^5 | 12 | 17 | 1024 → 256 | 13 |  |
| 3 | Thesis – C (with precomp. inverses & Fiat-Shamir) | 2^320 | 4096 | 32 | 5 | 15 | 2^17 | 2^5 | 12 | 17 | 1024 → 256 | 13 | Avg: 17.95 |
| 4 | Thesis – C (full optimisations) | 2^320 | 4096 | 32 | 5 | 15 | 2^17 | 2^5 | 12 | 17 | 1024 → 256 | 13 | Avg: 13.47 |
| 5 | Thesis – GPU (parallelised prover) | 2^320 | 4096 | 32 | 5 | 15 | 2^17 | 2^5 | 12 | 17 | 1024 → 256 | 14 | Avg: 9.36 |
| 6 | Thesis – GPU (with Merkle Tree) | 2^320 | 4096 | 32 | 5 | 15 | 2^17 | 2^5 | 12 | 17 | 1024 → 256 | 14 | Avg: 1.49 |

---

## Degree-Based Comparison

| S.No. | Degree Tested | Initial Length | Final Length | Time (Seconds) |
|-------|----------------|----------------|---------------|----------------|
| 5 | 512 | 2^14 | 2^5 | 0.543291 |
| 4 | 1024 | 2^15 | 2^5 | — |
| 3 | 2048 | 2^16 | 2^5 | 1.220363 |
| 1 | 4096 | 2^17 | 2^5 | 1.809525 |
| 2 | 8192 | 2^18 | 2^6 | 2.908838 |

---

## Sample Output Log (Excerpt)

```bash
Total non-zero basis elements: 200
Initial codeword length: 131072
Final codeword length after 12 rounds: 32
Number of rounds: 12
...
fri's prover took 0.677940 seconds to execute
Low degree test passed
Observed degree: 3
Should be no more than: 11
Verification successful!
fri took 1.408717 seconds to execute
Success! \o/
```
---

## Commands 
# Set shared library path (run every time)
export LD_LIBRARY_PATH=/u1/tjandhya/miniforge3/envs/sage/lib:$LD_LIBRARY_PATH

# Compile source files to object code
nvcc -dc ../src/file.cu -o file.o -Iinclude

# Link object files into executable
nvcc -o ../bin/fri merkle.o domain.o fft.o fiat-shamir.o params.o poly.o field.o commit-launch.o fri.o \
     -L/u1/tjandhya/miniforge3/envs/sage/lib -lcudart -lcrypto -Iinclude -cudart shared

# Run FRI prover/verifier
./fri

# Debug with GDB symbols enabled
nvcc -g -G -o ../bin/fri merkle.o domain.o fft.o fiat-shamir.o params.o poly.o field.o commit-launch.o fri.o \
     -L/u1/tjandhya/miniforge3/envs/sage/lib -lcudart -lcrypto -Iinclude -cudart shared

---
## Notes
	•	All timings are measured using time() from the <time.h> C library.
	•	Security levels L3 and L5 correspond to ~128-bit and ~192-bit post-quantum security respectively.
	•	The final round always compresses the codeword length to 32 or 64, depending on degree.
	•	All Merkle tree hashes are based on 1024-bit input → 256-bit output using Keccak or SHA-3.
