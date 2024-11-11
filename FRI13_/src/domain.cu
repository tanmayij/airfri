#include "../include/domain.cuh"
#include "../include/params.cuh"
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#define FIELD_WORDS preon.field_words
void initialize_domain(Domain *domain, size_t size, const uint64_t *basis, const uint64_t *shift) {
    assert(domain != NULL);

    domain->size = size;
    domain->basis_len = (size_t)log2(size);
    domain->basis = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));
    domain->shift = (uint64_t *)malloc(FIELD_WORDS * sizeof(uint64_t));

    memcpy((void *)domain->basis, (const void *)basis, FIELD_WORDS * sizeof(uint64_t));
    memcpy((void *)domain->shift, (const void *)shift, FIELD_WORDS * sizeof(uint64_t));
}

void free_domain(Domain *domain) {
    if (domain->basis) {
        free((void*)domain->basis);
        domain->basis = NULL;
    }
    if (domain->shift) {
        free((void*)domain->shift);
        domain->shift = NULL;
    }
}