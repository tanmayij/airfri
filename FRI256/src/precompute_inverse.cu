#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/field.cuh"  
#include "../include/params.cuh"

void calculate_inverses_to_file(const Parameters *params, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (file == NULL) {
        fprintf(stderr, "Error opening file: %s\n", filename);
        exit(EXIT_FAILURE);
    }

    // Loop over each FRI domain
    for (int d = 0; d < params->fri_num_reductions + 1; d++) {
        const Domain *domain = &params->fri_domains[d];
        size_t field_words = params->field_words; // e.g., 5 for 320-bit fields

        // Write the domain heading
        fprintf(file, "FRI Domain %d Inverses:\n", d);

        // Allocate memory for storing inverses of basis elements
        uint64_t **inverses = (uint64_t **)malloc(domain->basis_len * sizeof(uint64_t *));
        for (int i = 0; i < domain->basis_len; i++) {
            inverses[i] = (uint64_t *)malloc(field_words * sizeof(uint64_t));
        }

        // Calculate and write each inverse in the specified format
        for (int i = 0; i < domain->basis_len; i++) {
            const uint64_t *basis_element = &domain->basis[i * field_words];
            field_inv(inverses[i], basis_element, field_words);

            // Write each word of the inverse on a new line
            for (size_t w = 0; w < field_words; w++) {
                fprintf(file, "%016llx\n", (unsigned long long)inverses[i][w]);
            }
        }

        // Free allocated memory for inverses
        for (int i = 0; i < domain->basis_len; i++) {
            free(inverses[i]);
        }
        free(inverses);
    }

    fclose(file);
    printf("Inverses saved to %s\n", filename);
}

int main() {
    calculate_inverses_to_file(&preon, "fri_inverses_320_all.txt");
    return 0;
}