#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
// #include "../hash/hash.cuh"
//#include "../include/hash_host.cuh"
#include "../include/commit_kernel.cuh"
#include "../include/field.cuh"
//#include "../additive-fft/C++/Cantor/cantor_basis.hpp"
extern int field_words;
#include "../include/fri_utils.cuh"
using namespace std;
//hash
#include <stdio.h>
 #include <stdint.h>
 #include <string.h>
 
 #include "hash.cuh"
 
 #define SHA3_ASSERT( x )
 #define SHA3_TRACE( format, ...)
 #define SHA3_TRACE_BUF(format, buf, l)
 
 /* 
  * This flag is used to configure "pure" Keccak, as opposed to NIST SHA3.
  */
 #define SHA3_USE_KECCAK_FLAG 0x80000000
 #define SHA3_CW(x) ((x) & (~SHA3_USE_KECCAK_FLAG))
 
 #if defined(_MSC_VER)
 #define SHA3_CONST(x) x
 #else
 #define SHA3_CONST(x) x##L
 #endif
 
 #ifndef SHA3_ROTL64
 #define SHA3_ROTL64(x, y) \
     (((x) << (y)) | ((x) >> ((sizeof(uint64_t)*8) - (y))))
 #endif
 
 __device__ static const uint64_t keccakf_rndc[24] = {
     SHA3_CONST(0x0000000000000001UL), SHA3_CONST(0x0000000000008082UL),
     SHA3_CONST(0x800000000000808aUL), SHA3_CONST(0x8000000080008000UL),
     SHA3_CONST(0x000000000000808bUL), SHA3_CONST(0x0000000080000001UL),
     SHA3_CONST(0x8000000080008081UL), SHA3_CONST(0x8000000000008009UL),
     SHA3_CONST(0x000000000000008aUL), SHA3_CONST(0x0000000000000088UL),
     SHA3_CONST(0x0000000080008009UL), SHA3_CONST(0x000000008000000aUL),
     SHA3_CONST(0x000000008000808bUL), SHA3_CONST(0x800000000000008bUL),
     SHA3_CONST(0x8000000000008089UL), SHA3_CONST(0x8000000000008003UL),
     SHA3_CONST(0x8000000000008002UL), SHA3_CONST(0x8000000000000080UL),
     SHA3_CONST(0x000000000000800aUL), SHA3_CONST(0x800000008000000aUL),
     SHA3_CONST(0x8000000080008081UL), SHA3_CONST(0x8000000000008080UL),
     SHA3_CONST(0x0000000080000001UL), SHA3_CONST(0x8000000080008008UL)
 };
 
 __device__ static const unsigned keccakf_rotc[24] = {
     1, 3, 6, 10, 15, 21, 28, 36, 45, 55, 2, 14, 27, 41, 56, 8, 25, 43, 62,
     18, 39, 61, 20, 44
 };
 
 __device__ static const unsigned keccakf_piln[24] = {
     10, 7, 11, 17, 18, 3, 5, 16, 8, 21, 24, 4, 15, 23, 19, 13, 12, 2, 20,
     14, 22, 9, 6, 1
 };
 
 /* generally called after SHA3_KECCAK_SPONGE_WORDS-ctx->capacityWords words 
  * are XORed into the state s 
  */
  __device__ void keccakf(uint64_t s[25]) {
     int i, j, round;
     uint64_t t, bc[5];
     #define KECCAK_ROUNDS 24
 
     for (round = 0; round < KECCAK_ROUNDS; round++) {
         /* Theta */
         for (i = 0; i < 5; i++) {
             bc[i] = s[i] ^ s[i + 5] ^ s[i + 10] ^ s[i + 15] ^ s[i + 20];
         }
 
         for (i = 0; i < 5; i++) {
             t = bc[(i + 4) % 5] ^ SHA3_ROTL64(bc[(i + 1) % 5], 1);
             for (j = 0; j < 25; j += 5) {
                 s[j + i] ^= t;
             }
         }
 
         /* Rho Pi */
         t = s[1];
         for (i = 0; i < 24; i++) {
             j = keccakf_piln[i];
             bc[0] = s[j];
             s[j] = SHA3_ROTL64(t, keccakf_rotc[i]);
             t = bc[0];
         }
 
         /* Chi */
         for (j = 0; j < 25; j += 5) {
             for (i = 0; i < 5; i++) {
                 bc[i] = s[j + i];
             }
             for (i = 0; i < 5; i++) {
                 s[j + i] ^= (~bc[(i + 1) % 5]) & bc[(i + 2) % 5];
             }
         }
 
         /* Iota */
         s[0] ^= keccakf_rndc[round];
     }
 }
 
 /* *************************** Public Interface ************************ */
 
 /* For Init or Reset call these: */
 __device__ sha3_return_t sha3_Init(void *priv, unsigned bitSize) {
     sha3_context *ctx = (sha3_context *)priv;
     if (bitSize != 256 && bitSize != 384 && bitSize != 512) {
         return SHA3_RETURN_BAD_PARAMS;
     }
     memset(ctx, 0, sizeof(*ctx));  // Replace with cuda-compatible memset if needed
     ctx->capacityWords = 2 * bitSize / (8 * sizeof(uint64_t));
     return SHA3_RETURN_OK;
 }
 
 __device__ void sha3_Init256(void *priv) {
     sha3_Init(priv, 256);
 }
 
 __device__ void sha3_Init384(void *priv) {
     sha3_Init(priv, 384);
 }
 
 __device__ void sha3_Init512(void *priv) {
     sha3_Init(priv, 512);
 }
 
 __device__ enum SHA3_FLAGS sha3_SetFlags(void *priv, enum SHA3_FLAGS flags) {
     sha3_context *ctx = (sha3_context *)priv;
     flags = (enum SHA3_FLAGS)((uint32_t)flags & (uint32_t)SHA3_FLAGS_KECCAK);
     ctx->capacityWords |= (flags == SHA3_FLAGS_KECCAK ? SHA3_USE_KECCAK_FLAG : 0);
     return flags;
 }
 
 __device__ void sha3_Update(void *priv, void const *bufIn, size_t len) {
     sha3_context *ctx = (sha3_context *)priv;
 
     unsigned old_tail = (8 - ctx->byteIndex) & 7;
 
     size_t words;
     unsigned tail;
     size_t i;
 
     const uint8_t *buf = (const uint8_t *)bufIn;
 
     if (len < old_tail) {
         while (len--) {
             ctx->saved |= (uint64_t)(*(buf++)) << ((ctx->byteIndex++) * 8);
         }
         return;
     }
 
     if (old_tail) {
         len -= old_tail;
         while (old_tail--) {
             ctx->saved |= (uint64_t)(*(buf++)) << ((ctx->byteIndex++) * 8);
         }
         ctx->u.s[ctx->wordIndex] ^= ctx->saved;
         ctx->byteIndex = 0;
         ctx->saved = 0;
         if (++ctx->wordIndex == (SHA3_KECCAK_SPONGE_WORDS - SHA3_CW(ctx->capacityWords))) {
             keccakf(ctx->u.s);
             ctx->wordIndex = 0;
         }
     }
 
     words = len / sizeof(uint64_t);
     tail = len - words * sizeof(uint64_t);
 
     for (i = 0; i < words; i++, buf += sizeof(uint64_t)) {
         uint64_t t = (uint64_t)(buf[0]) |
                      ((uint64_t)(buf[1]) << 8 * 1) |
                      ((uint64_t)(buf[2]) << 8 * 2) |
                      ((uint64_t)(buf[3]) << 8 * 3) |
                      ((uint64_t)(buf[4]) << 8 * 4) |
                      ((uint64_t)(buf[5]) << 8 * 5) |
                      ((uint64_t)(buf[6]) << 8 * 6) |
                      ((uint64_t)(buf[7]) << 8 * 7);
         ctx->u.s[ctx->wordIndex] ^= t;
         if (++ctx->wordIndex == (SHA3_KECCAK_SPONGE_WORDS - SHA3_CW(ctx->capacityWords))) {
             keccakf(ctx->u.s);
             ctx->wordIndex = 0;
         }
     }
 
     while (tail--) {
         ctx->saved |= (uint64_t)(*(buf++)) << ((ctx->byteIndex++) * 8);
     }
 }
 
 __device__ void const *sha3_Finalize(void *priv) {
     sha3_context *ctx = (sha3_context *)priv;
     uint64_t t;
 
     if (ctx->capacityWords & SHA3_USE_KECCAK_FLAG) {
         t = (uint64_t)(((uint64_t)1) << (ctx->byteIndex * 8));
     } else {
         t = (uint64_t)(((uint64_t)(0x02 | (1 << 2))) << ((ctx->byteIndex) * 8));
     }
 
     ctx->u.s[ctx->wordIndex] ^= ctx->saved ^ t;
     ctx->u.s[SHA3_KECCAK_SPONGE_WORDS - SHA3_CW(ctx->capacityWords) - 1] ^= SHA3_CONST(0x8000000000000000UL);
     keccakf(ctx->u.s);
 
     return (ctx->u.sb);
 }
 
 __device__ sha3_return_t sha3_HashBuffer(unsigned bitSize, enum SHA3_FLAGS flags, const void *in, unsigned inBytes, void *out, unsigned outBytes) {
     sha3_return_t err;
     sha3_context c;
 
     err = sha3_Init(&c, bitSize);
     if (err != SHA3_RETURN_OK) {
         return err;
     }
     if (sha3_SetFlags(&c, flags) != flags) {
         return SHA3_RETURN_BAD_PARAMS;
     }
     sha3_Update(&c, in, inBytes);
     const void *h = sha3_Finalize(&c);
 
     if (outBytes > bitSize / 8) {
         outBytes = bitSize / 8;
     }
     memcpy(out, h, outBytes);
     return SHA3_RETURN_OK;
 }

__device__ void SHA3(uint8_t *hm, const uint8_t *msg, size_t msg_len, size_t bitSize)
{   
    SHA3s(hm, (const uint8_t *[]){msg}, (size_t[]){msg_len}, 1, bitSize);
}

__device__ void SHA3s(uint8_t *hm, const uint8_t **msgs, size_t *msgs_len, size_t items, size_t bitSize)
{
    sha3_context ctx;
    sha3_Init(&ctx, bitSize);
    for (size_t i = 0; i < items; ++i)
    {
        sha3_Update(&ctx, msgs[i], msgs_len[i]);
    }
    const void *result = sha3_Finalize(&ctx);
    memcpy(hm, result, bitSize / 8);
}
namespace cantor {
    constexpr uint64_t cantor_in_gf2to256[32][4] = {
{ 1ULL, 0ULL, 0ULL, 0ULL }, 
{ 4295033141ULL, 1ULL, 1ULL, 0ULL }, 
{ 13615767094782902071ULL, 16509767166769775667ULL, 17355192366725112762ULL, 17376248881838693815ULL }, 
{ 4644524726043978731ULL, 15693424397002034251ULL, 15710347196241526592ULL, 10880548238252614091ULL }, 
{ 640118059613985261ULL, 1313903195796726235ULL, 3396766875825798676ULL, 14192144114219446720ULL }, 
{ 10077849366122156232ULL, 7086092472166378657ULL, 15100903294859545647ULL, 3987845167958163125ULL }, 
{ 16270511978075044806ULL, 16164563811376147135ULL, 2625578537547985342ULL, 6072447074683886330ULL }, 
{ 5630018274122687930ULL, 8996055410665817774ULL, 3069990775405148075ULL, 12687162812642014724ULL }, 
{ 5754938134996921238ULL, 8937587289433511601ULL, 156331452736549519ULL, 10952228242332826081ULL }, 
{ 14899878026048491996ULL, 11881287613761720319ULL, 3157557464693671546ULL, 8490414338663093587ULL }, 
{ 6144948151697131770ULL, 12829366715966933238ULL, 5161617361959799949ULL, 16577793721757437839ULL }, 
{ 17335731826601981145ULL, 10374971554879084408ULL, 6560426406833767547ULL, 15376035997087156491ULL }, 
{ 18156021140703472090ULL, 4583267173085506923ULL, 8143265948793651304ULL, 13846076204877184192ULL }, 
{ 12368643584047939295ULL, 11819919576370109150ULL, 17178066221154598626ULL, 3537106026951133662ULL }, 
{ 17860052627717163084ULL, 1023859705452441292ULL, 4445661619295054034ULL, 16613301729520929541ULL }, 
{ 12567644037189678575ULL, 10989995650777694894ULL, 16650321402958599001ULL, 15394086554448135710ULL }, 
{ 3895293214564875702ULL, 6229113196406670081ULL, 3409072240043564264ULL, 9340437959674295503ULL }, 
{ 6003731047390282906ULL, 9123062914094717993ULL, 1995146645112670937ULL, 9677318193463286565ULL }, 
{ 11863262441345585988ULL, 16424244390609524150ULL, 12179255713357253792ULL, 8443460857916446437ULL }, 
{ 952441961121987866ULL, 17558130836699076420ULL, 13849621668659962270ULL, 5858758139358259385ULL }, 
{ 4842952427975513570ULL, 5536087762530766659ULL, 9066314258170351741ULL, 12689759026623861971ULL }, 
{ 11934042102125943131ULL, 5684306568414579395ULL, 388029345282030729ULL, 10845248503285669883ULL }, 
{ 3177192301120302685ULL, 15102579401182205913ULL, 12475081644470835988ULL, 14078415256170588310ULL }, 
{ 17763443710946410651ULL, 14346436287575459071ULL, 14335403489212834782ULL, 14214157082939204186ULL }, 
{ 5080490340599271256ULL, 12689939537085429179ULL, 1467244048409411438ULL, 14393358436686695531ULL }, 
{ 2614694687946290618ULL, 386260273944309489ULL, 3493066445529797887ULL, 14058644840294299294ULL }, 
{ 14488596929048962438ULL, 13785333696019006426ULL, 18177809138659256097ULL, 9536437080691907638ULL }, 
{ 9522776976361151227ULL, 12648439438049231636ULL, 9389018560787912221ULL, 9245451799559131985ULL }, 
{ 17882361982162706348ULL, 568081022572744249ULL, 13390524583313730048ULL, 14385592888705335992ULL }, 
{ 4104017265333605612ULL, 15082008458465896366ULL, 16686233106678210090ULL, 3681558827891653969ULL }, 
{ 7929606058812646144ULL, 3003101923230410791ULL, 13498978919646001267ULL, 11922546952913825908ULL }, 
{ 18231145936351310928ULL, 11618316067220172972ULL, 7652517275783664224ULL, 2804331428984460444ULL }, 
};
}

#define HASH_WORDS 4
#define FIELD_WORDS 4
#define CONCAT_WORDS (FIELD_WORDS + HASH_WORDS)


inline void cantor_Lk(uint64_t* out, const uint64_t basis[][4], int basis_len, size_t idx) {
    for(int w = 0; w < 4; w++) out[w] = 0;
    for(int b = 0; b < basis_len; b++) {
        if((idx >> b) & 1) {
            for(int w = 0; w < 4; w++) {
                out[w] ^= basis[b][w];
            }
        }
    }
}

__device__ void print_field_kernel(const char *label, const uint64_t *field, int field_words) {
    printf("%s: ", label);
    for (int i = 0; i < field_words; i++) {
        printf("%016lx ", field[i]);
    }
    printf("\n");
}

__host__ void print_field_host(const char *label, const uint64_t *field, int field_words) {
    printf("%s: ", label);
    for (int i = 0; i < field_words; i++) {
        printf("%016lx ", field[i]);
    }
    printf("\n");
}

__host__ void initialize_file(const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return;
    }
    fclose(file); 
}

__host__ void write_to_file(const char *filename, const uint64_t *data, int field_words, int total_indices) {
    FILE *file = fopen(filename, "a");
    if (file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        return;
    }

    for (int index = 0; index < total_indices; index++) {
        fprintf(file, "Index %d: ", index);
        for (int i = 0; i < field_words; i++) {
            fprintf(file, "%llx ", (unsigned long long)data[index + i]);
        }
        fprintf(file, "\n");
    }
    fflush(file);  
    fclose(file);
}

__global__ void commit_kernel(
    uint64_t *device_codeword, uint64_t *device_codeword_nxt,
    uint64_t *device_alpha,
    uint64_t *device_temp1, uint64_t *device_temp2,
    uint64_t *device_temp3, uint64_t *device_temp4,
    uint64_t *device_temp5,
    uint64_t *device_alpha_offset, 
    int N
) {
    size_t I = blockIdx.x * blockDim.x + threadIdx.x;

    if (I >= N / 2) return;

    int idx1 = 2 * I * FIELD_WORDS;
    int idx2 = (2 * I + 1) * FIELD_WORDS;
    int idx3 = I * FIELD_WORDS;

    // Cantor basis: f_{k+1}(L_{k+1}[i]) = (f_k(L_k[2i]) - f_k(L_k[2i + 1])) * (alpha - L_k[2i]) + f_k(L_k[2i])
    // Let L_k[2i] = x_even, L_k[2i+1] = x_odd
    // temp1 = codeword[2j] - codeword[2j+1]
    field_sub(&device_temp1[idx3], &device_codeword[idx1], &device_codeword[idx2], FIELD_WORDS);
    // temp2 = alpha - x_even (L_k[2i])
    field_sub(&device_temp2[idx3], device_alpha,  &device_alpha_offset[idx3], FIELD_WORDS);
    // temp3 = temp1 * temp2
    field_mul(&device_temp3[idx3], &device_temp1[idx3], &device_temp2[idx3], FIELD_WORDS);
    // temp4 = codeword[2j] + temp3
    field_add(&device_temp4[idx3], &device_codeword[idx1], &device_temp3[idx3], FIELD_WORDS);
    memcpy(&device_codeword_nxt[idx3], &device_temp4[idx3], FIELD_WORDS * sizeof(uint64_t));
}


__global__ void compute_tree_layers(
    uint64_t *device_codeword_nxt, 
    uint64_t *device_layer_hashes, 
    uint64_t *device_tree_layer,
    uint64_t *device_tree_layer_nxt, 
    uint64_t *device_combined_sibling_codewords, 
    uint64_t *device_concat_codeword_to_hash, 
    uint64_t *device_digest, 
    int N
) {   
    size_t I = blockIdx.x * blockDim.x + threadIdx.x;

    if (I >= N / 2) return;

    if (I < N / 2 && N == 1048576) {
        int idx1 = 2 * I * FIELD_WORDS;
        int idx2 = (2 * I + 1) * FIELD_WORDS;
        int idx3 = I * (2 * FIELD_WORDS);  // for codeword_nxt element
        int idx5 = I * HASH_WORDS;   // for hash
        int idx4 = I * (FIELD_WORDS + HASH_WORDS);

        for (int j = 0; j < FIELD_WORDS; j++) {
            device_combined_sibling_codewords[idx1 + j] = device_tree_layer[idx1 + j];
            device_combined_sibling_codewords[idx1 + FIELD_WORDS + j] = device_tree_layer[idx2 + j];
        }

        SHA3((uint8_t *)&device_digest[idx5], (uint8_t *)&device_combined_sibling_codewords[idx3], 2 * FIELD_WORDS * sizeof(uint64_t), 256);

        for (int j = 0; j < FIELD_WORDS; j++) {
            device_concat_codeword_to_hash[idx4 + j] = device_codeword_nxt[idx1 / 2 + j]; //fix
        }
        for (int j = 0; j < HASH_WORDS; j++) {
            device_concat_codeword_to_hash[idx4 + FIELD_WORDS + j] = device_digest[idx5 + j];
        }
        for (int j = 0; j < (FIELD_WORDS + HASH_WORDS); j++) {
            device_tree_layer_nxt[idx4 + j] = device_concat_codeword_to_hash[idx4 + j];
        }
    }

    if (I < N / 2 && N < 1048576 && N >= 64) {
        int idx1 = 2 * I * (FIELD_WORDS + HASH_WORDS);
        int idx2 = ((2 * I) + 1) * (FIELD_WORDS + HASH_WORDS);
        int idx3 = I * FIELD_WORDS;
        int idx5 = I * HASH_WORDS;
        int idx4 = I * (FIELD_WORDS + HASH_WORDS);
        int idx6 = I * 2 * (FIELD_WORDS + HASH_WORDS);

        for (int j = 0; j < (FIELD_WORDS + HASH_WORDS); j++) {
            device_combined_sibling_codewords[idx6 + j] = device_tree_layer[idx1 + j];
            device_combined_sibling_codewords[idx6 + (FIELD_WORDS + HASH_WORDS) + j] = device_tree_layer[idx2 + j];
        }

        SHA3((uint8_t *)&device_digest[idx5], (uint8_t *)&device_combined_sibling_codewords[idx6], 2 * (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t), 256);

    
        for (int j = 0; j < FIELD_WORDS; j++) {
            device_concat_codeword_to_hash[idx4 + j] = device_codeword_nxt[idx3 + j];  //fixed??
        }
        for (int j = 0; j < HASH_WORDS; j++) {
            device_concat_codeword_to_hash[idx4 + FIELD_WORDS + j] = device_digest[idx5 + j];
        }

        for (int j = 0; j < (FIELD_WORDS + HASH_WORDS); j++) {
            device_tree_layer_nxt[idx4 + j] = device_concat_codeword_to_hash[idx4 + j];
        }
    }
}
__global__ void merkle_kernel(
    uint64_t *device_layer_hashes, 
    uint64_t *device_merkle_root, 
    uint64_t *device_tree_layer,
    uint64_t *device_tree_layer_nxt,
    uint64_t *device_combined_sibling_codewords,
    uint64_t *device_digest,
    uint64_t *device_combined_sibling_hashes,
    int N
) {
    int I = blockIdx.x * blockDim.x + threadIdx.x;
    int idx1, idx2, idx3, idx4, idx5;

    if (I >= N / 2) return;

    if (N == 32) {  // Only if N is 32, we use codewords
        idx1 = 2 * I * (FIELD_WORDS + HASH_WORDS);
        idx2 = (2 * I + 1) * (FIELD_WORDS + HASH_WORDS);
        idx3 = I * HASH_WORDS;
        idx4 = I * (2 * (FIELD_WORDS + HASH_WORDS));
        idx5 = I * (FIELD_WORDS + HASH_WORDS);

        for (int j = 0; j < (FIELD_WORDS + HASH_WORDS); j++) {
            device_combined_sibling_codewords[idx4 + j] = device_tree_layer[idx1 + j];
            device_combined_sibling_codewords[idx4 + (FIELD_WORDS + HASH_WORDS) + j] = device_tree_layer[idx2 + j];
        }

        SHA3((uint8_t *)&device_digest[idx3], (uint8_t *)&device_combined_sibling_codewords[idx4], 2 * (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t), 256);

        for (int j = 0; j < HASH_WORDS; j++) {
            device_tree_layer_nxt[idx3 + j] = device_digest[idx3 + j];
        }
    }

    if (N < 32) {  
        idx1 = 2 * I * HASH_WORDS;
        idx2 = (2 * I + 1) * HASH_WORDS;
        idx3 = I * HASH_WORDS;
        idx4 = I * (2 * HASH_WORDS);
        
        for (int j = 0; j < HASH_WORDS; j++) {
            device_combined_sibling_hashes[idx4 + j] = device_tree_layer[idx1 + j];
            device_combined_sibling_hashes[idx4 + HASH_WORDS + j] = device_tree_layer[idx2 + j];
        }
        printf("combined_sibling_hashes before SHA3 in merkle kernel: ");
        for (int j = 0; j < 2 * HASH_WORDS; j++) {
            printf("%016lx ", device_combined_sibling_hashes[j]);
        }
        printf("\n");

        SHA3((uint8_t *)&device_digest[idx3], (uint8_t *)&device_combined_sibling_hashes[idx4], 2 * HASH_WORDS * sizeof(uint64_t), 256);

        for (int j = 0; j < HASH_WORDS; j++) {
            device_tree_layer_nxt[idx3 + j] = device_digest[idx3 + j];
        }
    }
}


__global__ void compute_merkle_root_kernel(
    uint64_t *device_tree_layer,    // Input: layer with two sibling hashes
    uint64_t *device_merkle_root    // Output: the Merkle root
) {
    if (threadIdx.x == 0) {
        int idx1 = 0 * HASH_WORDS; // First sibling hash
        int idx2 = 1 * HASH_WORDS; // Second sibling hash

        uint64_t combined_sibling_hashes[2 * HASH_WORDS];

        for (int j = 0; j < HASH_WORDS; j++) {
            combined_sibling_hashes[j] = device_tree_layer[idx1 + j];
            combined_sibling_hashes[HASH_WORDS + j] = device_tree_layer[idx2 + j];
        }

        printf("combined_sibling_hashes before SHA3 in commit kernel: ");
        for (int j = 0; j < 2 * HASH_WORDS; j++) {
            printf("%016lx ", combined_sibling_hashes[j]);
        }
        printf("\n");

        SHA3((uint8_t *)device_merkle_root, (uint8_t *)combined_sibling_hashes, 
             2 * HASH_WORDS * sizeof(uint64_t), 256);

        printf("Computed Merkle Root inside kernel: ");
        for (int j = 0; j < HASH_WORDS; j++) {
            printf("%016lx ", device_merkle_root[j]);
        }
        printf("\n");
    }
}
// __global__ void compute_merkle_root_kernel(
//     uint64_t *device_tree_layer,    // Input: layer with two sibling hashes
//     uint64_t *device_merkle_root   // Output: the Merkle root
// ) {
//     if (threadIdx.x == 0) {
//         // Indices for the two sibling hashes
//         int idx1 = 0 * HASH_WORDS;
//         int idx2 = 1 * HASH_WORDS;

//         // Combined buffer for the two sibling hashes
//         uint64_t combined_sibling_hashes[2 * HASH_WORDS];

//         // Copy the sibling hashes into the combined buffer
//         memcpy(&combined_sibling_hashes[0], &device_tree_layer[idx1], HASH_WORDS * sizeof(uint64_t));
//         memcpy(&combined_sibling_hashes[HASH_WORDS], &device_tree_layer[idx2], HASH_WORDS * sizeof(uint64_t));
//         printf("combined_sibling_hashes before SHA3 in commit kernel: ");
//         for (int j = 0; j < 2 * HASH_WORDS; j++) {
//             printf("%016lx ", combined_sibling_hashes[j]);
//         }
//         printf("\n");
//         // Compute the Merkle root using SHA3
//         SHA3(
//             (uint8_t *)device_merkle_root,                  // Destination: Merkle root
//             (uint8_t *)combined_sibling_hashes,            // Source: Combined sibling hashes
//             2 * HASH_WORDS * sizeof(uint64_t),             // Input size: 2 sibling hashes
//             256                                            // Output size: 256 bits
//         );
//         printf("Computed Merkle Root inside kernel: ");
//         for (int j = 0; j < HASH_WORDS; j++) {
//             printf("%016lx ", device_merkle_root[j]);
//         }
//         printf("\n");
//     }
// }

void commit_launch(
    uint64_t **codeword, uint64_t **codeword_nxt, 
    uint64_t *alpha, int N, uint64_t *root, uint64_t **tree_layer, uint64_t **tree_layer_nxt, uint64_t ***tree, int last_round, bool is_last_round
) {
    printf("Starting commit_launch\n");
    printf("N = %d, FIELD_WORDS = %d\n", N, FIELD_WORDS);
    int basis_len = (int)log2(N);
    printf("basis len: %d\n", basis_len);
    printf("last round: %d\n", last_round);
    if (N == 1048576) {
        initialize_file("temp1.txt");
        initialize_file("temp2.txt");
        initialize_file("temp3.txt");
        initialize_file("temp4.txt");
        initialize_file("temp5.txt");
    }

    uint64_t *flattened_codeword = (uint64_t *)malloc(N * FIELD_WORDS * sizeof(uint64_t));
    if (flattened_codeword == NULL) {
        fprintf(stderr, "Error: malloc failed for flattened_codeword\n");
        return;
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < FIELD_WORDS; ++j) {
            int index = i * FIELD_WORDS + j;
            flattened_codeword[index] = codeword[i][j];
        }
    }

    printf("First few flattened codeword values:\n");
    for (int i = 0; i < 10; ++i) {
        printf("%016lx ", flattened_codeword[i]);
    }
    printf("\n");

    uint64_t *flattened_tree_layer;
    if(N == 1048576) {
    flattened_tree_layer = (uint64_t *)malloc(N * FIELD_WORDS * sizeof(uint64_t));
    if (flattened_tree_layer == NULL) {
        fprintf(stderr, "Error: malloc failed for flattened_codeword\n");
        return;
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < FIELD_WORDS; ++j) {
            int index = i * FIELD_WORDS + j;
            flattened_tree_layer[index] = tree_layer[i][j];
        }
    }
    printf("First few flattened tree_layer values:\n");
    for (int i = 0; i < 16; ++i) {
        printf("%016lx ", flattened_tree_layer[i]);
    }
    printf("\n");

    } else {
    flattened_tree_layer = (uint64_t *)malloc(N * CONCAT_WORDS * sizeof(uint64_t));
    if (flattened_tree_layer == NULL) {
        fprintf(stderr, "Error: malloc failed for flattened_codeword\n");
        return;
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < CONCAT_WORDS; ++j) {
            int index = i * CONCAT_WORDS + j;
            flattened_tree_layer[index] = tree_layer[i][j];
        }
    }
    printf("First few flattened tree_layer values:\n");
    for (int i = 0; i < 16; ++i) {
        printf("%016lx ", flattened_tree_layer[i]);
    }
    printf("\n");
}


    uint64_t *flattened_codeword_nxt = (uint64_t *)malloc((N / 2) * FIELD_WORDS * sizeof(uint64_t));
    if (flattened_codeword_nxt == NULL) {
        fprintf(stderr, "Error: malloc failed for flattened_codeword_nxt\n");
        free(flattened_codeword);
        return;
    }

    for (int i = 0; i < N / 2; ++i) {
        codeword_nxt[i] = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
        if (codeword_nxt[i] == NULL) {
            fprintf(stderr, "Error: malloc failed for codeword_nxt[%d]\n", i);
            for (int j = 0; j < i; ++j) {
                free(codeword_nxt[j]);
            }
            free(flattened_codeword);
            free(flattened_codeword_nxt);
            return;
        }
    }

    uint64_t *flattened_tree_layer_nxt = (uint64_t *)malloc((N / 2) * CONCAT_WORDS * sizeof(uint64_t));
    if (flattened_tree_layer_nxt == NULL) {
        fprintf(stderr, "Error: malloc failed for flattened_tree_layer_nxt\n");
        free(flattened_codeword);
        return;
    }

    for (int i = 0; i < N / 2; ++i) {
        tree_layer_nxt[i] = (uint64_t *)malloc(CONCAT_WORDS * sizeof(uint64_t));
        if (tree_layer_nxt[i] == NULL) {
            fprintf(stderr, "Error: malloc failed for tree_layer_nxt[%d]\n", i);
            for (int j = 0; j < i; ++j) {
                free(tree_layer_nxt[j]);
            }
            free(flattened_tree_layer_nxt);
            return;
        }
    }


    int field_size = N * FIELD_WORDS * sizeof(uint64_t);

    uint64_t *device_codeword, *device_codeword_nxt, *device_alpha;
    uint64_t *device_temp1, *device_temp2, *device_temp3, *device_temp4, *device_temp5;
    uint64_t *device_layer_hashes, *device_merkle_root, *device_tree_layer, *device_tree_layer_nxt, *device_combined_sibling_codewords, *device_combined_sibling_hashes, *device_concat_codeword_to_hash, *device_digest;
    //uint64_t *flattened_alpha_offset = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp1 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp2 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp3 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp4 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t *flattened_temp5 = (uint64_t *)malloc(N/2 * FIELD_WORDS * sizeof(uint64_t));
    uint64_t* flattened_alpha_offset = (uint64_t*)malloc((N/2) * FIELD_WORDS * sizeof(uint64_t));
    int cantor_basis_len = basis_len; // basis_len for this round
    for(int i = 0; i < N/2; i++) {
        cantor_Lk(&flattened_alpha_offset[i * FIELD_WORDS], cantor::cantor_in_gf2to256, cantor_basis_len, 2 * i);
    }
    uint64_t* device_alpha_offset;
    cudaMalloc((void**)&device_alpha_offset, (N/2) * FIELD_WORDS * sizeof(uint64_t));
    
    flattened_temp1[N/2 * FIELD_WORDS] = {0}, flattened_temp2[N/2 * FIELD_WORDS] = {0}, flattened_temp3[N/2 * FIELD_WORDS] = {0}, flattened_temp4[N/2 * FIELD_WORDS] = {0}, flattened_temp5[N/2 * FIELD_WORDS] = {0};
    flattened_alpha_offset[N/2 * FIELD_WORDS] = {0};
    cudaMalloc((void**)&device_codeword, field_size);
    cudaMalloc((void**)&device_codeword_nxt, (N/2) * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void**)&device_alpha, FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void **)&device_temp1, N/2 * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void **)&device_temp2, N/2 * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void **)&device_temp3, N/2 * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void **)&device_temp4, N/2 * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void **)&device_temp5, N/2 * FIELD_WORDS * sizeof(uint64_t));
    cudaMalloc((void**)&device_layer_hashes, N * HASH_WORDS * sizeof(uint64_t));
    if(N == 1048576){
    cudaMalloc((void**)&device_tree_layer, field_size);
    cudaMalloc((void**)&device_combined_sibling_codewords, (N/2) * 2 * FIELD_WORDS * sizeof(uint64_t));
    } else 
    {
        cudaMalloc((void**)&device_tree_layer, N * CONCAT_WORDS * sizeof(uint64_t));
        cudaMalloc((void**)&device_combined_sibling_codewords, (N/2) * 2 * ( FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
    }
    cudaMalloc((void**)&device_tree_layer_nxt, (N/2) * CONCAT_WORDS * sizeof(uint64_t));
    cudaMalloc((void**)&device_merkle_root, HASH_WORDS * sizeof(uint64_t));
    cudaMalloc((void**)&device_concat_codeword_to_hash, (N/2) * (FIELD_WORDS + HASH_WORDS) * sizeof(uint64_t));
    cudaMalloc((void**)&device_combined_sibling_hashes, (N/2) * 2 * HASH_WORDS * sizeof(uint64_t));
    cudaMalloc((void**)&device_digest, (N/2) * HASH_WORDS * sizeof(uint64_t));

    cudaMemcpy(device_codeword, flattened_codeword, field_size, cudaMemcpyHostToDevice);
    cudaMemcpy(device_alpha, alpha, FIELD_WORDS * sizeof(uint64_t), cudaMemcpyHostToDevice);
    cudaMemcpy(device_alpha_offset, flattened_alpha_offset, (N/2) * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyHostToDevice);
    // cudaMemcpy(device_offset, offset, sizeof(uint64_t), cudaMemcpyHostToDevice);
    // cudaMemcpy(device_denominator_inv, &denominator_inv, sizeof(uint64_t), cudaMemcpyHostToDevice);
    // cudaMemcpy(device_eval_basis, flattened_eval_basis, basis_len * sizeof(uint64_t), cudaMemcpyHostToDevice);
    if(N == 1048576) {
    cudaMemcpy(device_tree_layer, flattened_tree_layer,  N * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyHostToDevice);
    } else 
    {
        cudaMemcpy(device_tree_layer, flattened_tree_layer,  N * CONCAT_WORDS * sizeof(uint64_t), cudaMemcpyHostToDevice);
    }
    int threads_per_block = 1;
    int num_blocks = (N / 2 + threads_per_block - 1) / threads_per_block;
    commit_kernel<<<num_blocks * 2, threads_per_block>>>(
        device_codeword, device_codeword_nxt, device_alpha,
        device_temp1, device_temp2, device_temp3, device_temp4, device_temp5, device_alpha_offset, N
    );
    cudaDeviceSynchronize();

    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Kernel launch failed: %s\n", cudaGetErrorString(err));
        free(flattened_codeword);
        free(flattened_codeword_nxt);
        // free(flattened_tree_layer);
        //free(flattened_eval_basis);
        return;
    }
    
    compute_tree_layers<<<num_blocks * 2, threads_per_block>>> (
        device_codeword_nxt, device_layer_hashes, device_tree_layer, device_tree_layer_nxt, device_combined_sibling_codewords, device_concat_codeword_to_hash, device_digest, N
    );
    cudaDeviceSynchronize();

    err = cudaGetLastError();
    if (err != cudaSuccess) {
        printf("Kernel launch failed: %s\n", cudaGetErrorString(err));
        free(flattened_codeword);
        free(flattened_codeword_nxt);
        // free(flattened_tree_layer);
        //free(flattened_eval_basis);
        return;
    }

    cudaMemcpy(flattened_codeword_nxt, device_codeword_nxt, (N / 2) * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp1, device_temp1, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp2, device_temp2, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp3, device_temp3, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp4, device_temp4, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_temp5, device_temp5, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    // cudaMemcpy(flattened_alpha_offset, device_alpha_offset, N/2 * FIELD_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
    cudaMemcpy(flattened_tree_layer_nxt, device_tree_layer_nxt, (N / 2) * CONCAT_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);

    

    write_to_file("temp1.txt", flattened_temp1, FIELD_WORDS, N/2);
    write_to_file("temp2.txt", flattened_temp2, FIELD_WORDS, N/2);
    write_to_file("temp3.txt", flattened_temp3, FIELD_WORDS, N/2);
    write_to_file("temp4.txt", flattened_temp4, FIELD_WORDS, N/2);
    write_to_file("temp5.txt", flattened_temp5, FIELD_WORDS, N/2);
    // write_to_file("alpha_offset.txt", flattened_alpha_offset, FIELD_WORDS, N/2);
    
    for (int i = 0; i < N / 2; ++i) {
        for (int j = 0; j < FIELD_WORDS; ++j) {
            int index = i * FIELD_WORDS + j;
            if (index >= (N / 2) * FIELD_WORDS) {
                printf("Out-of-bounds access at index: %d\n", index);
                break;
            }
            codeword_nxt[i][j] = flattened_codeword_nxt[index];
        }
    }

    printf("First few codeword_nxt values:\n");
    for (int i = 0; i < 10; i++) {
        printf("%016lx ", flattened_codeword_nxt[i]);
    }
    printf("\n");
    cudaFree(device_codeword);

    for (int i = 0; i < N / 2; ++i) {
        for (int j = 0; j < CONCAT_WORDS; ++j) {
            int index = i * CONCAT_WORDS + j;
            if (index >= (N / 2) * CONCAT_WORDS) {
                printf("Out-of-bounds access at index: %d\n", index);
                break;
            }
            tree_layer_nxt[i][j] = flattened_tree_layer_nxt[index];
        }
    }

    printf("First few tree_layer_nxt values:\n");
    for (int i = 0; i < 10; i++) {
        printf("%016lx ", flattened_tree_layer_nxt[i]);
    }
    printf("\n");
    if(!is_last_round){
    free(flattened_tree_layer_nxt);
    }
    cudaFree(device_tree_layer);
    cudaFree(device_tree_layer_nxt);

    if (is_last_round) { 
        
        int tree_idx = last_round;  //start with tree[15] which will hold 32 elements
        //pring last_round value
        cout << "Last round value in commit_launch: " << last_round << endl;
        int next_N = N / 2; //initialize to 32 for the next layer size

        //first transfer the tree_layer_nxt elements to tree[15]
        tree[tree_idx] = (uint64_t **)malloc((next_N) * sizeof(uint64_t *));
        for (int i = 0; i < next_N; i++) {
            tree[tree_idx][i] = (uint64_t *)malloc(CONCAT_WORDS * sizeof(uint64_t));
            for (int j = 0; j < CONCAT_WORDS; j++) {
                tree[tree_idx][i][j] = flattened_tree_layer_nxt[i * CONCAT_WORDS + j];
            }
        }
        printf("\n=== DEBUG: flattened_tree_layer_nxt After Storing in tree[%d] ===\n", tree_idx); //to check if the elements have codewords in them
        for (int i = 0; i < next_N; i++) {
            printf("Index %d: ", i);
            for (int j = 0; j < CONCAT_WORDS; j++) {
                printf("%016lx ", flattened_tree_layer_nxt[i * CONCAT_WORDS + j]);
            }
            printf("\n");
        }
        /*steps:
        1. tree[14] is filled, flattened_tree_layer_nxt has 32 elements now
        2. allocate memory to device variables tree and tree_nxt (host to device)
        3.copy device_tree_layer = flattened_tree_layer_nxt
        4.realloc memory for flattened_tree_layer and flattened_tree_layer_nxt
        5. then launch merkle_kernel  */
    
        //copy the tree_layer_nxt to tree_layer
        cudaMalloc((void **)&device_tree_layer, (next_N) * CONCAT_WORDS * sizeof(uint64_t)); //32
        cudaMalloc((void **)&device_tree_layer_nxt, (next_N / 2) * HASH_WORDS * sizeof(uint64_t)); //will be computed in the next few lines //16
        cudaMemcpy(device_tree_layer, flattened_tree_layer_nxt, (next_N) * CONCAT_WORDS * sizeof(uint64_t), cudaMemcpyHostToDevice); //32
    
        //now, flattened_tree_layer will be 32 and flattened_tree_layer_nxt will be 16
        flattened_tree_layer = (uint64_t *)realloc(flattened_tree_layer, (next_N) * CONCAT_WORDS * sizeof(uint64_t)); //32
        flattened_tree_layer_nxt = (uint64_t *)realloc(flattened_tree_layer_nxt, (next_N / 2) * HASH_WORDS * sizeof(uint64_t)); //16
       
        tree_idx++;  //move to the next tree layer index, which is 1

        // // Step 3: transfer 
        // cudaMemcpy(device_tree_layer, device_tree_layer_nxt, (next_N / 2) * CONCAT_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToDevice);
        // merkle_kernel<<<2, 32>>>(
        //     device_layer_hashes, 
        //     device_merkle_root, 
        //     device_tree_layer, 
        //     device_tree_layer_nxt, 
        //     device_combined_sibling_codewords, 
        //     device_digest, 
        //     device_combined_sibling_hashes, 
        //     next_N
        // );
        // cudaDeviceSynchronize();
        // //transfer device_tree_layer here 
        // cudaMemcpy(device_tree_layer, device_tree_layer_nxt, (next_N / 2) * HASH_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToDevice);

        //step 3: Loop over remaining layers, updating tree[layer] with each iteration
        while (next_N > 2) {
            //print current next_N 
            cout << "Current next_N in loop: " << next_N << endl;
            int tpb = min(32, next_N / 2);
            if(tpb == 0) {tpb = 1;}
            int nb = (next_N + tpb - 1) / tpb;

            // merkle_kernel for each layer (computes next layer hashes)
            merkle_kernel<<<nb, tpb>>>(device_layer_hashes, device_merkle_root, device_tree_layer, device_tree_layer_nxt, device_combined_sibling_codewords, device_digest, device_combined_sibling_hashes, next_N);
            cudaDeviceSynchronize();
            
            //cudaMalloc((void **)&device_tree_layer_nxt, (next_N / 2) * HASH_WORDS * sizeof(uint64_t));
            cudaMemcpy(flattened_tree_layer_nxt, device_tree_layer_nxt, (next_N / 2) * HASH_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
            printf("First few tree_layer_nxt values in the looooop:\n");
            for (int i = 0; i < HASH_WORDS * 4; i++) {
                printf("%016lx ", flattened_tree_layer_nxt[i]);
            }
            printf("\n");
            //unflatten and store in tree[tree_idx]
            printf("Populated tree[%d] with %d elements\n", tree_idx, next_N/2);
            tree[tree_idx] = (uint64_t **)malloc((next_N / 2) * sizeof(uint64_t *));
            for (int i = 0; i < next_N / 2; i++) {
                tree[tree_idx][i] = (uint64_t *)malloc(HASH_WORDS * sizeof(uint64_t));
                for (int j = 0; j < HASH_WORDS; j++) {
                    tree[tree_idx][i][j] = flattened_tree_layer_nxt[i * HASH_WORDS + j];
                }
            }
            cudaFree(device_combined_sibling_codewords);
            cudaMalloc((void**)&device_combined_sibling_codewords, (next_N / 2) * 2 * HASH_WORDS * sizeof(uint64_t)); //pick just hashes after merkle_launch launches the first time. 
            uint64_t *temp = device_tree_layer; 
            device_tree_layer = device_tree_layer_nxt;
            // device_tree_layer_nxt = temp;
            tree_idx++;
            //print tree_idx
            cout << "Current tree_idx in loop: " << tree_idx << endl;
            next_N /= 2;
            //print next_N
            cout << "Updated next_N in loop: " << next_N << endl;
            flattened_tree_layer_nxt = (uint64_t *)realloc(flattened_tree_layer_nxt, (next_N) * HASH_WORDS * sizeof(uint64_t)); //fixed from CONCAT_WORDS to HASH_WORDS
            //step 4: Update device_tree_layer with the contents of device_tree_layer_nxt
            //cudaMemcpy(device_tree_layer, device_tree_layer_nxt, (next_N / 2) * HASH_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToDevice);


        }// close while
        //do not delete the below code - only comment it
        int number = 0;
        for (int layer = 15; layer <= 19; layer++) {
            int elements = 1 << (20 - layer);  // Number of elements in this layer
            printf("\n=== DEBUG: tree[%d] Elements ===\n", layer);
            if (layer == 15) { number = 8;}
            else { number = 4;}
            for (int i = 0; i < elements; i++) {
                printf("Index %d: ", i);
                for (int j = 0; j < number; j++) {
                    printf("%016lx ", tree[layer][i][j]);  // Print each hash
                }
                printf("\n");
            }
        }
        printf("\n=== DEBUG: tree[18] Elements (Last Layer Before Merkle Root) ===\n");
        for (int i = 0; i < 2; i++) {  // 2 elements in layer 18
            printf("Index %d: ", i);
            for (int j = 0; j < 4; j++) {
                printf("%016lx ", tree[18][i][j]);  // Print each hash
            }
            printf("\n");
        }

        compute_merkle_root_kernel<<<1, 1>>>(
            device_tree_layer,  // The last computed layer (contains 2 hashes)
            device_merkle_root      // Output: The final Merkle root
        );
        cudaDeviceSynchronize();
        
        cudaMemcpy(root, device_merkle_root, HASH_WORDS * sizeof(uint64_t), cudaMemcpyDeviceToHost);
        
        printf("Computed Merkle Root: ");
        for (int i = 0; i < HASH_WORDS; i++) {
            printf("%016lx ", root[i]);
        }
        printf("\n");
        
    cudaFree(device_tree_layer);
    cudaFree(device_tree_layer_nxt);
    cudaFree(device_codeword_nxt);
    cudaFree(device_alpha);
    cudaFree(device_merkle_root);
    cudaFree(device_layer_hashes);
    }
    free(flattened_codeword);
    free(flattened_codeword_nxt);


    printf("Memory freed and commit_launch completed.\n");
}
