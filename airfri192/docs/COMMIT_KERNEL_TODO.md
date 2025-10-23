# Commit Kernel Implementation Notes

## What I've Created

### Files:
1. **`src/commit_kernel.cu`** - Main CUDA implementation
2. **`include/commit_kernel.hpp`** - C++ header for host code

## Current Status

### ✅ Completed:
- Basic kernel structure for FRI folding
- Device-compatible field operations (add, sub, mul placeholders)
- `i_th_element_in_span` for Cantor basis
- Merkle tree layer computation structure
- Host launch function skeleton
- **libff::gf256 ↔ uint64_t[4] conversion** (`include/field_conversion.hpp`)
- **SHA3 device functions integrated** (using `hash.cu`)
- **Unified field operations** (using proven `field.cu` from FRI192)

### ⚠️ TODO - Critical Items:

#### 1. **✅ libff::gf256 to uint64_t Conversion (COMPLETED)**
**Solution**: Created `include/field_conversion.hpp` with:
- `field_to_uint64()`: libff::gf256 → uint64_t[4]
- `uint64_to_field()`: uint64_t[4] → libff::gf256
- Uses libff's `convert_field_element_to_bit_vector()` internally

#### 2. **✅ Proper GF(2^256) Multiplication (COMPLETED)**
**Solution**: Using proven `field.cu` implementation
- All field operations unified: `field_add`, `field_sub`, `field_mul`, `field_inv`
- `__host__ __device__` compatible
- Same code works on CPU and GPU
- Tested implementation from FRI192

**Options:**
a) Port libff's gf256::operator* to CUDA
b) Copy your CUDA fri.cu field_mul implementation
c) Use lookup tables (if small enough)

**Recommendation**: Copy from `FRI192/include/field.cuh::field_mul()`

#### 3. **✅ SHA3 Device Function (COMPLETED)**
**Solution**: Using existing `hash.cu` with SHA3 device functions
- Updated `commit_kernel.cu` to call `SHA3()` directly
- Removed placeholder XOR code
b) Use existing CUDA SHA3 library
c) Call host SHA3 (slow, not recommended)

**Where**: Line ~174 in `compute_tree_layer` kernel

#### 4. **Basis Inverse Computation**
**Current**: `flat_basis_inv[FIELD_WORDS] = {0}` (line ~242)

**Need**: Precompute `inverse(basis[0])` on host
- Use libff's `.inverse()` method
- Only needed once per round

#### 5. **Debug File Writing**
**I commented out** all file writes to match your request

**To enable debugging:**
Uncomment lines ~130-140 (in commit_kernel)
```cpp
//if(I == 0 && blockIdx.x == 0) {
//    print_field_device("result", result);
//}
```

## Architecture Notes

### Layer Structure (15 FRI rounds, 2^20 → 2^5):

```
Layer  | Size   | Content         | Kernel
-------|--------|-----------------|--------
0      | 2^20   | field only      | N/A (input)
1      | 2^19   | field || hash   | commit + tree
2      | 2^18   | field || hash   | commit + tree
...
15     | 2^5    | field || hash   | commit + tree
16     | 2^4    | hash only       | merkle only
17     | 2^3    | hash only       | merkle only
18     | 2^2    | hash only       | merkle only
19     | 2^1    | hash only       | merkle only
20     | 2^0    | hash only (root)| merkle only
```

### FRI Folding Formula:
```
next[j] = current[2j] + alpha * (current[2j+1] - current[2j]) / (alpha - domain[2j])
```

**Simplification with Cantor basis:**
- Domain has special structure
- Can use precomputed inverse of basis[0]

## Integration with prove.cpp

```cpp
//in commit_host():
for(int r = 0; r < 15; r++) {
    //hash alpha for this round
    //...
    
    //launch GPU kernel
    commit_launch(codeword, codeword_nxt, alpha, offset, basis,
                 N, r, merkle_tree_layers);
    
    //update for next round
    codeword = codeword_nxt;
    N = N / 2;
}
```

## Next Steps

1. **Implement libff ↔ uint64_t conversion**
   - Check libff/gf256.hpp for internal structure
   - May need `.data` member or `.as_ulong_long()` method

2. **Port field multiplication**
   - Copy from FRI192/include/field.cuh
   - Add `__device__` tag
   - Test with known inputs

3. **Port SHA3**
   - Copy from FRI192/include/hash.cuh  
   - Adapt for device code
   - May need separate .cuh file

4. **Test incrementally**
   - Test field ops first (add, mul)
   - Test kernel with small N (32)
   - Verify against CPU version
   - Scale up to 2^20

5. **Optimize**
   - Shared memory for basis access
   - Coalesced memory access
   - Bank conflict reduction

## Debugging Tips

```cpp
//uncomment in kernel to see values:
if(I == 0) {
    print_field_device("alpha", device_alpha);
    print_field_device("result", result);
}

//on host:
printf("first few codeword_nxt: ");
for(int i = 0; i < 10; i++) {
    printf("%016lx ", flat_codeword_nxt[i]);
}
```

## Notes on Your Requirements

✅ **Flattening**: Using `uint64_t* flat_array` with `[i * FIELD_WORDS + w]` indexing
✅ **15 rounds**: Updated from 12 to 15 (2^20 → 2^5)  
✅ **Layer structure**: 0=field, 1-15=field||hash, 16-20=hash
✅ **No denominator_inv param**: Removed, using basis_inv instead
⚠️ **libff functions**: Need to port or extract data manually
✅ **Small comments**: Using `//comment` style
✅ **Debugging**: Commented out, easy to uncomment

## Complexity Warning

**Using different field operations WILL cause debugging headaches!**

**Recommendation**: 
1. Create `field_device.cuh` with exact copies of your CUDA field ops
2. Add `__device__` tags
3. Test thoroughly
4. Then integrate with libff types only at boundaries (host code)

This way you have:
- **Host**: Clean C++ with libff types
- **Device**: Raw CUDA with uint64_t[4]
- **Boundary**: Conversion functions

