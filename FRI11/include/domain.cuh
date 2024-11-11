#ifndef _DOMAIN_H__
#define _DOMAIN_H__

#include <stddef.h>
#include <stdint.h>

typedef struct Domain {
    size_t size;
    size_t basis_len; // is log2(size)
    const uint64_t *basis;
    const uint64_t *shift;
} Domain;

void initialize_domain(Domain *domain, size_t size, const uint64_t *basis, const uint64_t *shift);
void free_domain(Domain *domain);

#endif // _DOMAIN_H__